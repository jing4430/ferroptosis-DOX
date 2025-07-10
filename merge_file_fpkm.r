# 加载所需的包
library(data.table)

# 定义文件路径列表
file_paths <- c("GSE76314.txt","GSE106297.txt","GSE157282.txt","GSE247345.txt")  # 替换为实际文件路径

# 初始化一个空 data.table，用于存储合并后的数据
combined_data <- NULL

# 遍历文件路径，分块读取并合并每个文件的数据
for (file in file_paths) {
  # 打开文件连接
  con <- file(file, open = "r")
  
  # 读取文件的第一行（列名）
  header <- readLines(con, n = 1)
  header <- unlist(strsplit(header, "\t"))
  
  # 检查是否包含 gene_symbol 列
  if (!"gene_symbol" %in% header) {
    stop(paste("文件", file, "中缺少 gene_symbol 列"))
  }
  
  # 初始化一个空 data.table，用于存储当前文件的数据
  file_data <- NULL
  
  # 分块读取文件内容
  chunk_size <- 10000  # 每次读取 10000 行
  while (TRUE) {
    chunk <- readLines(con, n = chunk_size)
    if (length(chunk) == 0) break  # 如果读取完毕，退出循环
    
    # 将块数据转换为 data.table
    chunk_data <- fread(text = chunk, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
    setnames(chunk_data, header)  # 设置列名
    
    # 在同一个文件内，对相同的 gene_symbol 进行合并，其他列的值取平均数
    chunk_data <- chunk_data[, lapply(.SD, mean, na.rm = TRUE), by = gene_symbol]
    
    # 合并块数据
    if (is.null(file_data)) {
      file_data <- chunk_data
    } else {
      file_data <- rbindlist(list(file_data, chunk_data))[, lapply(.SD, mean, na.rm = TRUE), by = gene_symbol]
    }
    
    # 释放临时变量以节省内存
    rm(chunk_data)
    gc()  # 强制垃圾回收
  }
  
  # 关闭文件连接
  close(con)
  
  # 合并当前文件的数据到 combined_data
  if (is.null(combined_data)) {
    combined_data <- file_data
  } else {
    combined_data <- merge(combined_data, file_data, by = "gene_symbol", all = TRUE)
  }
  
  # 释放临时变量以节省内存
  rm(file_data)
  gc()  # 强制垃圾回收
}

# 删除包含空值的行
combined_data <- na.omit(combined_data)

# 打印合并后的数据
print(combined_data)

# 将合并后的数据写入新的 .txt 文件
fwrite(combined_data, file = "merged_output.txt", sep = "\t", na = "NA", quote = FALSE)

# 提示合并完成
cat("文件合并完成，结果已保存到 merged_output.txt\n")
