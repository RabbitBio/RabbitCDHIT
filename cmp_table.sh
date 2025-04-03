# chmod +x cmp_table.sh
# ./cmp_table.sh

cmp -s master_word_table.info worker_word_table.info
if [ $? -eq 0 ]; then
    echo "文件完全一致"
else
    echo "文件不一致"
fi
