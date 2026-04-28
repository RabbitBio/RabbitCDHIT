#!/bin/bash
# 单机一步运行 RabbitCD-HIT
# 用法: bash run_single.sh -i <input.fa> -o <output> [选项]
#
# 示例:
#   bash run_single.sh -i DB.fa -o result -c 0.9 -T 64

set -e

# ── 环境变量 ──────────────────────────────────────────────────
export OMP_PLACES=cores
export OMP_PROC_BIND=close
export I_MPI_PIN=1
export I_MPI_DEBUG=5

# ── 默认参数 ──────────────────────────────────────────────────
INPUT=""
OUTPUT=""
TOTAL_CORES=$(nproc)                                                      # 整机逻辑核心数
NUMA_NODES=$(LC_ALL=C lscpu | awk '/^NUMA node\(s\)/{print $3}')          # NUMA 节点数
PREPROCESS_T=$((NUMA_NODES > 0 ? TOTAL_CORES / NUMA_NODES : 32))         # preprocess -T：单 NUMA 核心数，探测失败默认 32
# MPI_T 不在此处设置，必须从 preprocess 生成的 info.json 中读取（见下方）
CLUSTER_THD=0.9     # 相似度阈值 -c
MEMORY_CONTROL=""   # memory_control: "" = 自动检测, "1" = 强制开启, "0" = 强制关闭
EXTRA_PREPROCESS="" # 传给 cdhit-preprocess 的额外参数
EXTRA_MPI=""        # 传给 cdhit-mpi 的额外参数

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
TMP_DIR="${SCRIPT_DIR}/tmp_runs"
PRE_OUT="${SCRIPT_DIR}/preprocess_output"

# ── 解析参数 ──────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
    case "$1" in
        -i)   INPUT="$2";      shift 2 ;;
        -o)   OUTPUT="$2";     shift 2 ;;
        -T)   PREPROCESS_T="$2"; shift 2 ;;  # 仅控制 preprocess 线程数，mpi -T 固定从 info.json 读取
        -c)   CLUSTER_THD="$2"; shift 2 ;;
        -tmp) TMP_DIR="$2";    shift 2 ;;
        -pre_out) PRE_OUT="$2"; shift 2 ;;
        -memory_control) MEMORY_CONTROL="$2"; shift 2 ;;
        *)
            # 其余参数根据前缀分发给对应阶段
            # 以 -stealing / -g / -fo / -G 等开头的给 mpi
            case "$1" in
                -stealing|-g|-fo|-G|-b|-n|-t|-s|-S|-aL|-AL|-aS|-AS|-A|-uL|-uS|-U|-d|-p|-sc|-sf|-bak)
                    EXTRA_MPI="$EXTRA_MPI $1 $2"; shift 2 ;;
                *)
                    echo "未知参数: $1" >&2; exit 1 ;;
            esac
            ;;
    esac
done

# ── 参数检查 ──────────────────────────────────────────────────
if [[ -z "$INPUT" || -z "$OUTPUT" ]]; then
    echo "用法: bash $0 -i <input.fa> -o <output> [-T 线程数] [-c 相似度] [其他选项]"
    echo ""
    echo "  -i        输入 FASTA 文件（支持 .gz）"
    echo "  -o        输出文件前缀"
    echo "  -T        preprocess 线程数，默认自动取单 NUMA 核心数（cdhit-mpi 的 -T 固定从 info.json 读取，不受此参数影响）"
    echo "  -c        相似度阈值，默认 0.9"
    echo "  -tmp      临时文件目录，默认 ./tmp_runs"
    echo "  -pre_out  预处理输出目录，默认 ./preprocess_output"
    echo "  -memory_control <0|1>  内存控制模式："
    echo "              不传    自动检测（默认）：内存不足时自动开启"
    echo "              1       强制开启：按需从磁盘加载序列，降低峰值内存"
    echo "              0       强制关闭：禁用自动检测，所有序列常驻内存"
    echo "            注意：-memory_control 1 与 -stealing 不可同时使用"
    echo "  其他 cdhit-mpi 参数（-g, -fo, -stealing, -G, -b, -n, -t, -s, -S ...）直接附加即可"
    exit 1
fi

INPUT="$(realpath "$INPUT")"
TMP_DIR="$(realpath "$TMP_DIR")"
PRE_OUT="$(realpath "$PRE_OUT")"

# ── Step 1: 预处理 ────────────────────────────────────────────
echo "========================================"
echo "Step 1: 预处理"
echo "  输入:       $INPUT"
echo "  整机核数(-NT): $TOTAL_CORES"
echo "  单NUMA核数(-T): $PREPROCESS_T  (NUMA节点数: $NUMA_NODES)"
echo "  临时目录:   $TMP_DIR"
echo "  预处理输出: $PRE_OUT"
echo "========================================"

"${SCRIPT_DIR}/cdhit-preprocess" \
    -i "$INPUT" \
    -N 1 \
    -NT "$TOTAL_CORES" \
    -T  "$PREPROCESS_T" \
    -c  "$CLUSTER_THD" \
    -tmp    "$TMP_DIR" \
    -pre_out "$PRE_OUT"

# ── 从 info.json 读取 MPI 参数 ────────────────────────────────
INFO_JSON="${PRE_OUT}/info.json"
if [[ ! -f "$INFO_JSON" ]]; then
    echo "错误: 找不到 ${INFO_JSON}，预处理可能失败" >&2
    exit 1
fi

# cdhit-mpi 的 -np 和 -T 必须与 info.json 严格一致，不可手动覆盖
TOTAL_MPI=$(python3 -c "import json; d=json.load(open('${INFO_JSON}')); print(d['info']['total_mpi_num'])")
MPI_T=$(python3 -c "import json; d=json.load(open('${INFO_JSON}')); print(d['info']['threads_per_rank'])")

echo ""
echo "========================================"
echo "Step 2: 聚类"
echo "  mpirun -np $TOTAL_MPI  (from info.json: total_mpi_num)"
echo "  cdhit-mpi -T $MPI_T    (from info.json: threads_per_rank)"
if [[ -n "$MEMORY_CONTROL" ]]; then
    echo "  -memory_control $MEMORY_CONTROL"
else
    echo "  -memory_control      (auto-detect)"
fi
echo "  输出前缀: $OUTPUT"
echo "========================================"

# 若用户指定了 memory_control，拼入 MPI 参数
MC_ARG=""
if [[ -n "$MEMORY_CONTROL" ]]; then
    MC_ARG="-memory_control $MEMORY_CONTROL"
fi

mpirun -np "$TOTAL_MPI" \
    "${SCRIPT_DIR}/cdhit-mpi" \
    -i  "$PRE_OUT" \
    -o  "$OUTPUT"  \
    -T  "$MPI_T" \
    $MC_ARG \
    $EXTRA_MPI

echo ""
echo "========================================"
echo "完成！输出文件:"
echo "  ${OUTPUT}      (代表序列)"
echo "  ${OUTPUT}.clstr（聚类信息）"
echo "========================================"
