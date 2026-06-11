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
LOAD_ALL=""         # load_all: "" = 双缓冲流式（默认）, "1" = 全量加载到内存
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
        -load_all) LOAD_ALL="$2"; shift 2 ;;
        *)
            # 其余参数根据前缀分发给对应阶段
            # 以 -stealing / -g / -fo / -G 等开头的给 mpi
            case "$1" in
                -stealing|-g|-fo|-G|-b|-n|-t|-s|-S|-aL|-AL|-aS|-AS|-A|-uL|-uS|-U|-d|-p|-sc|-sf|-bak)
                    EXTRA_MPI="$EXTRA_MPI $1 $2"; shift 2 ;;
                *)
                    echo "Unknown argument: $1" >&2; exit 1 ;;
            esac
            ;;
    esac
done

# ── 参数检查 ──────────────────────────────────────────────────
if [[ -z "$INPUT" || -z "$OUTPUT" ]]; then
    echo "Usage: bash $0 -i <input.fa> -o <output> [-T threads] [-c identity] [other options]"
    echo ""
    echo "  -i        Input FASTA file (.gz supported)"
    echo "  -o        Output file prefix"
    echo "  -T        Number of threads for preprocessing. By default, it is automatically set to the number of cores per NUMA node. The -T value for cdhit-mpi is read from info.json and is not affected by this option."
    echo "  -c        Sequence identity threshold. Default: 0.9"
    echo "  -tmp      Temporary directory. Default: ./tmp_runs"
    echo "  -pre_out  Preprocessing output directory. Default: ./preprocess_output"
    echo "  -load_all <0|1>  Sequence loading mode:"
    echo "              not set / 0  Double-buffered streaming (default): sequences are read from disk"
    echo "                           per chunk and freed after use; saves memory for large datasets"
    echo "              1            Load all: keep all sequences resident in RAM"
    echo "            Note: -load_all 0 (streaming) cannot be combined with -stealing"
    echo "  Other cdhit-mpi options (-g, -fo, -stealing, -G, -b, -n, -t, -s, -S ...) can be appended directly"
    exit 1
fi

INPUT="$(realpath "$INPUT")"
TMP_DIR="$(realpath "$TMP_DIR")"
PRE_OUT="$(realpath "$PRE_OUT")"

# ── Step 1: 预处理 ────────────────────────────────────────────
echo "========================================"
echo "Step 1: Preprocessing"
echo "  Input:                 $INPUT"
echo "  Total cores (-NT):      $TOTAL_CORES"
echo "  Cores per NUMA (-T):    $PREPROCESS_T  (NUMA nodes: $NUMA_NODES)"
echo "  Temporary directory:    $TMP_DIR"
echo "  Preprocessing output:   $PRE_OUT"
echo "========================================"

"${SCRIPT_DIR}/cdhit-preprocess" \
    -i "$INPUT" \
    -N 1 \
    -NT "$TOTAL_CORES" \
    -T  "$PREPROCESS_T" \
    -tmp    "$TMP_DIR" \
    -pre_out "$PRE_OUT"

# ── 从 info.json 读取 MPI 参数 ────────────────────────────────
INFO_JSON="${PRE_OUT}/info.json"
if [[ ! -f "$INFO_JSON" ]]; then
    echo "Error: ${INFO_JSON} was not found. Preprocessing may have failed." >&2
    exit 1
fi

# cdhit-mpi 的 -np 和 -T 必须与 info.json 严格一致，不可手动覆盖
TOTAL_MPI=$(python3 -c "import json; d=json.load(open('${INFO_JSON}')); print(d['info']['total_mpi_num'])")
MPI_T=$(python3 -c "import json; d=json.load(open('${INFO_JSON}')); print(d['info']['threads_per_rank'])")

echo ""
echo "========================================"
echo "Step 2: Clustering"
echo "  mpirun -np $TOTAL_MPI  (from info.json: total_mpi_num)"
echo "  cdhit-mpi -T $MPI_T    (from info.json: threads_per_rank)"
if [[ -n "$LOAD_ALL" ]]; then
    echo "  -load_all $LOAD_ALL"
else
    echo "  -load_all            (default: 0, double-buffered streaming)"
fi
echo "  Output prefix: $OUTPUT"
echo "========================================"

# 若用户显式指定了 load_all，拼入 MPI 参数
MC_ARG=""
if [[ -n "$LOAD_ALL" ]]; then
    MC_ARG="-load_all $LOAD_ALL"
fi

mpirun -np "$TOTAL_MPI" \
    "${SCRIPT_DIR}/cdhit-mpi" \
    -i  "$PRE_OUT" \
    -o  "$OUTPUT"  \
    -T  "$MPI_T" \
    -c  "$CLUSTER_THD" \
    $MC_ARG \
    $EXTRA_MPI

echo ""
echo "========================================"
echo "Done. Output files:"
echo "  ${OUTPUT}        (representative sequences)"
echo "  ${OUTPUT}.clstr  (clustering information)"
echo "========================================"