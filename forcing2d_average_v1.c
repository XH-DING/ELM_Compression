/*
 * 并行读写NetCDF文件的示例：forcing2d_average
 * 功能：根据指定的变量类型、年份和月份，找到对应的NetCDF文件
 * M个进程按time维度分割读取一个文件，对自己的部分求平均
 * 然后进程间取平均，最后得到[y,x]大小的数据
 * 写入时按y维度分割，每个进程负责一部分
 */

 #include <stdio.h>
 #include <stdlib.h>
 #include <string.h>
 #include <dirent.h>
 #include <mpi.h>
 #include <pnetcdf.h>
 #include <unistd.h>  /* 用于getopt */
 
 /* 错误处理宏 */
 #define CHECK_ERR(err) { \
     if (err != NC_NOERR) { \
         printf("Error at line %d: %s\n", __LINE__, ncmpi_strerror(err)); \
         MPI_Abort(MPI_COMM_WORLD, -1); \
         return 1; \
     } \
 }
 
 /* 最大文件路径长度 */
 #define MAX_PATH_LEN 1024
 /* 最大文件数量 */
 #define MAX_FILES 1000
 /* 最大变量类型数量 */
 #define MAX_VAR_TYPES 100
 
 // int xlen_nc_type(nc_type xtype, int *size)
 // {
 //     switch(xtype) {
 //         case NC_BYTE:
 //         case NC_CHAR:
 //         case NC_UBYTE:  *size = 1; return NC_NOERR;
 //         case NC_SHORT:
 //         case NC_USHORT: *size = 2; return NC_NOERR;
 //         case NC_INT:
 //         case NC_UINT:
 //         case NC_FLOAT:  *size = 4; return NC_NOERR;
 //         case NC_DOUBLE:
 //         case NC_INT64:
 //         case NC_UINT64: *size = 8; return NC_NOERR;
 //         default: return NC_EBADTYPE;
 //     }
 // }
 
 
 /* 查找匹配的文件 */
 int find_matching_files(const char *input_dir, const char **var_types, int num_var_types, 
                         int year, int month, char **file_list, int *num_files) {
     DIR *dir;
     struct dirent *entry;
     char pattern[256];
     int found = 0;
     
     dir = opendir(input_dir);
     if (dir == NULL) {
         printf("Error: Could not open directory %s\n", input_dir);
         return -1;
     }
     
     *num_files = 0;
     
     /* 对每个变量类型进行匹配 */
     for (int i = 0; i < num_var_types; i++) {
         /* 构建匹配模式: clmforc.Daymet4.1km.VARTYPE.YYYY-MM.nc */
         sprintf(pattern, "clmforc.Daymet4.1km.%s.%04d-%02d.nc", var_types[i], year, month);
         
         /* 重置目录读取位置 */
         rewinddir(dir);
         
         /* 遍历目录寻找匹配文件 */
         while ((entry = readdir(dir)) != NULL) {
             if (strcmp(entry->d_name, pattern) == 0) {
                 /* 构建完整文件路径 */
                 sprintf(file_list[*num_files], "%s/%s", input_dir, entry->d_name);
                 (*num_files)++;
                 found = 1;
                 break; /* 找到一个就跳出，避免重复添加 */
             }
         }
     }
     
     closedir(dir);
     return found ? 0 : -1;
 }
 
 /* 解析逗号分隔的变量列表 */
 int parse_variable_list(const char *var_string, char **var_types, int *num_var_types) {
     char *token;
     char *string_copy = strdup(var_string);
     char *save_ptr = NULL;
     int count = 0;
     
     token = strtok_r(string_copy, ",", &save_ptr);
     while (token != NULL && count < MAX_VAR_TYPES) {
         /* 去除前后空格 */
         char *start = token;
         while (*start == ' ') start++;
         
         char *end = start + strlen(start) - 1;
         while (end > start && *end == ' ') end--;
         *(end + 1) = '\0';
         
         var_types[count] = strdup(start);
         count++;
         
         token = strtok_r(NULL, ",", &save_ptr);
     }
     
     *num_var_types = count;
     free(string_copy);
     return 0;
 }
 
 /* 显示使用帮助 */
 void show_usage(const char *program_name) {
     printf("Usage: %s [options]\n", program_name);
     printf("Options:\n");
     printf("  -i <input_path>  指定输入路径\n");
     printf("  -o <output_path> 指定输出文件路径(文件名将自动生成为forcing2d_average_YYYY_MM.nc)\n");
     printf("  -y <year>        指定年份\n");
     printf("  -m <month>       指定月份\n");
     printf("  -v <variables>   指定变量列表，以逗号分隔\n");
     printf("  -h               显示帮助信息\n");
     printf("Example:\n");
     printf("  mpiexec -n 14 %s -i /input/path -o output/path -y 2014 -m 12 -v FLDS,FSDS,WIND\n", program_name);
 }
 
 int main(int argc, char **argv) {
     int ret, i, j;
     // int nc_size;
     int ncid_in, ncid_out, varid_in, *varid_out;
     int *dimids_in, *dimids_out, ndims;
     MPI_Offset *dim_sizes_in, *dim_sizes_out;
     int global_rank, global_size;
     int file_group, proc_in_group, num_groups, procs_per_group;
     char **input_files;
     char output_path[MAX_PATH_LEN] = "";
     char output_file[MAX_PATH_LEN];
     float *buffer = NULL;           // 输入缓冲区
     float *local_avg = NULL;        // 局部平均值
     float *global_avg = NULL;       // 全局平均值
     MPI_Comm file_comm;
     MPI_Info info;
     char **var_types;               // 变量类型数组
     
     /* 添加计时变量 */
     double start_time, end_time, process_time;
     double read_start, read_end, read_time;
     double compute_start, compute_end, compute_time;
     double write_start_time, write_end_time, write_time;
     double total_read_time, total_compute_time, total_write_time, total_time;
 
     /* 参数相关变量 */
     char input_dir[MAX_PATH_LEN] = "";
     int num_var_types = 0;
     int year = -1;
     int month = -1;
     int num_files = 0;
     char var_string[MAX_PATH_LEN] = "";
     int opt;
     
 
     /* 初始化MPI */
     MPI_Init(&argc, &argv);
     MPI_Comm_rank(MPI_COMM_WORLD, &global_rank);
     MPI_Comm_size(MPI_COMM_WORLD, &global_size);
     
     /* 开始计时 - 整个程序开始 */
     start_time = MPI_Wtime();
 
     /* 解析命令行参数 */
     while ((opt = getopt(argc, argv, "i:o:y:m:v:h")) != -1) {
         switch (opt) {
             case 'i':
                 strcpy(input_dir, optarg);
                 break;
             case 'o':
                 strcpy(output_path, optarg);
                 break;
             case 'y':
                 year = atoi(optarg);
                 break;
             case 'm':
                 month = atoi(optarg);
                 break;
             case 'v':
                 strcpy(var_string, optarg);
                 break;
             case 'h':
                 if (global_rank == 0) {
                     show_usage(argv[0]);
                 }
                 MPI_Finalize();
                 return 0;
             default:
                 if (global_rank == 0) {
                     fprintf(stderr, "Unknown option: %c\n", opt);
                     show_usage(argv[0]);
                 }
                 MPI_Finalize();
                 return 1;
         }
     }
     
     /* 检查必要参数 */
     if (input_dir[0] == '\0' || output_path[0] == '\0' || 
         year < 0 || month < 1 || month > 12 || var_string[0] == '\0') {
         if (global_rank == 0) {
             fprintf(stderr, "Error: Missing required parameters\n");
             show_usage(argv[0]);
         }
         MPI_Finalize();
         return 1;
     }
     
     /* 构建输出文件名 */
     /* 检查输出路径是否以斜杠结尾 */
     if (output_path[strlen(output_path) - 1] != '/') {
         sprintf(output_file, "%s/forcing2d_average_%04d_%02d.nc", output_path, year, month);
     } else {
         sprintf(output_file, "%sforcing2d_average_%04d_%02d.nc", output_path, year, month);
     }
 
     /* 解析变量列表 */
     var_types = (char **)malloc(MAX_VAR_TYPES * sizeof(char *));
     parse_variable_list(var_string, var_types, &num_var_types);
     
     if (num_var_types == 0) {
         if (global_rank == 0) {
             fprintf(stderr, "Error: No valid variables specified\n");
         }
         MPI_Finalize();
         return 1;
     }
     
     if (global_rank == 0) {
         printf("输入目录: %s\n", input_dir);
         printf("输出文件: %s\n", output_file);
         printf("年份: %d\n", year);
         printf("月份: %d\n", month);
         printf("变量列表 (%d个): ", num_var_types);
         for (i = 0; i < num_var_types; i++) {
             printf("%s%s", var_types[i], (i < num_var_types - 1) ? ", " : "\n");
         }
     }
     
     /* 分配文件列表内存 */
     input_files = (char **)malloc(MAX_FILES * sizeof(char *));
     for (i = 0; i < MAX_FILES; i++) {
         input_files[i] = (char *)malloc(MAX_PATH_LEN * sizeof(char));
     }
     
     /* 查找匹配的文件 */
     if (global_rank == 0) {
         printf("开始查找匹配的文件...\n");
         ret = find_matching_files(input_dir, (const char **)var_types, num_var_types, year, month, input_files, &num_files);
         if (ret != 0 || num_files == 0) {
             printf("Error: No matching files found\n");
             MPI_Abort(MPI_COMM_WORLD, -1);
             return 1;
         }
         printf("找到 %d 个匹配的文件\n", num_files);
         for (i = 0; i < num_files; i++) {
             printf("文件 %d: %s\n", i, input_files[i]);
         }
     }
     
     /* 将文件数广播给所有进程 */
     MPI_Bcast(&num_files, 1, MPI_INT, 0, MPI_COMM_WORLD);
     
     /* 将文件名广播给所有进程 */
     for (i = 0; i < num_files; i++) {
         MPI_Bcast(input_files[i], MAX_PATH_LEN, MPI_CHAR, 0, MPI_COMM_WORLD);
     }
     
     
     /* 检查进程数是否是文件数的倍数 */
     if (global_size % num_files != 0) {
         if (global_rank == 0) {
             printf("Warning: Number of processes (%d) is not a multiple of the number of files (%d)\n", global_size, num_files);
             printf("Some processes may remain idle\n");
         }
     }
     /* 检查进程数是否大于文件数 */
     if (global_size / num_files == 0) {
         if (global_rank == 0) {
             printf("Error: Number of processes (%d) is less than the number of files (%d)\n", global_size, num_files);
             MPI_Finalize();
             return 0;
         }
     }
 
     /* 设置组数为文件数 */
     num_groups = num_files;
 
     /* 计算每组的进程数 */
     procs_per_group = global_size / num_groups;
     
     /* 确定当前进程属于哪个文件组和组内的序号 */
     file_group = global_rank / procs_per_group;
     proc_in_group = global_rank % procs_per_group;
     
     /* 如果进程所属组超出了文件数，则该进程不参与计算 */
     if (file_group >= num_files) {
         MPI_Finalize();
         return 0;
     }
     
     /* 创建当前进程组的通信域 */
     MPI_Comm_split(MPI_COMM_WORLD, file_group, proc_in_group, &file_comm);
     
     /* 创建MPI信息对象 */
     MPI_Info_create(&info);
     MPI_Info_set(info, "nc_chunk_default_filter", "sz");
     MPI_Info_set(info, "nc_chunking", "enable");
     
     /* 第一阶段：读取输入文件 */
     if (global_rank == 0) {
         printf("开始读取输入文件...\n");
     }
     
     /* 开始读取计时 */
     read_start = MPI_Wtime();
 
     /* 打开当前组对应的输入文件 */
     ret = ncmpi_open(file_comm, input_files[file_group], NC_NOWRITE, info, &ncid_in);
     CHECK_ERR(ret);
     
    //  printf("1111");
     /* 获取变量ID - 使用文件名中的变量类型名作为变量名 */
     char var_name[NC_MAX_NAME+1];
     strcpy(var_name, var_types[file_group]);
     ret = ncmpi_inq_varid(ncid_in, var_name, &varid_in);
     CHECK_ERR(ret);
    //  printf("2222");
     /* 获取变量维度数 */
     ret = ncmpi_inq_varndims(ncid_in, varid_in, &ndims);
     CHECK_ERR(ret);
    //  printf("3333");
     if (ndims != 3) {
         printf("Error: Expected 3 dimensions (time, y, x) but found %d dimensions\n", ndims);
         MPI_Abort(MPI_COMM_WORLD, -1);
         return 1;
     }
    //  printf("4444");
     /* 分配维度ID和大小数组 */
     dimids_in = (int *)malloc(ndims * sizeof(int));
     dim_sizes_in = (MPI_Offset *)malloc(ndims * sizeof(MPI_Offset));
     
     /* 获取变量的维度ID */
     ret = ncmpi_inq_vardimid(ncid_in, varid_in, dimids_in);
     CHECK_ERR(ret);
    //  printf("5555");
     /* 获取每个维度的大小和名称 */
     char **dim_names = (char **)malloc(ndims * sizeof(char *));
    //  for (i = 0; i < ndims; i++) {
    //     printf("dimids_in[%d]: %d\n", i, dimids_in[i]);
    //  }
     for (i = 0; i < ndims; i++) {
        //  dim_names[i] = (char *)malloc((NC_MAX_NAME+1) * sizeof(char));
        //  printf("6666");
        //  ret = ncmpi_inq_dimname(ncid_in, dimids_in[i], dim_names[i]);
        //  CHECK_ERR(ret);
        //  printf("7777");
         ret = ncmpi_inq_dimlen(ncid_in, dimids_in[i], &dim_sizes_in[i]);
        //  printf("dim_size[%d]: %d\n", i, dim_sizes_in[i]);
         CHECK_ERR(ret);
     }
     dim_names[0] = "time";
     dim_names[1] = "y";
     dim_names[2] = "x";
     /* 计算空间维度大小（y * x）*/
     MPI_Offset spatial_size = dim_sizes_in[1] * dim_sizes_in[2]; // y * x
     MPI_Offset time_steps = dim_sizes_in[0]; // time dimension size
     
     /* 设置读取起始位置和计数 - 在time维度上分割 */
     MPI_Offset time_chunk = time_steps / procs_per_group;
     MPI_Offset time_remainder = time_steps % procs_per_group;
     MPI_Offset my_time_count = (proc_in_group < time_remainder) ? time_chunk + 1 : time_chunk;
     MPI_Offset my_time_start = (proc_in_group < time_remainder) ? proc_in_group * (time_chunk + 1) : proc_in_group * time_chunk + time_remainder;
     
     /* 分配读取起始位置和计数数组 */
     MPI_Offset start[3], count[3];
     
     /* 在time维度上分割 */
     start[0] = my_time_start;
     count[0] = my_time_count;
     /* 读取完整的y和x维度 */
     start[1] = 0;
     count[1] = dim_sizes_in[1];
     start[2] = 0;
     count[2] = dim_sizes_in[2];
     
     /* 获取变量的数据类型 */
     // nc_type var_type;
     // ret = ncmpi_inq_vartype(ncid_in, varid_in, &var_type);
     // CHECK_ERR(ret);
     // ret = xlen_nc_type(var_type, nc_size)
     // CHECK_ERR(ret);
    // printf("8888");
     /* 分配内存用于读取数据 */
     MPI_Offset local_elements = my_time_count * spatial_size;
     buffer = (float *)malloc(local_elements * sizeof(float));
     if (buffer == NULL) {
         printf("Error: Memory allocation failed for buffer\n");
         MPI_Abort(MPI_COMM_WORLD, -1);
         return 1;
     }
    //  printf("9999");
     /* 读取数据 */
     ret = ncmpi_get_vara_float_all(ncid_in, varid_in, start, count, buffer);
     CHECK_ERR(ret);
     
     /* 关闭输入文件 */
     ret = ncmpi_close(ncid_in);
     CHECK_ERR(ret);
    //  printf("AAAA");
     /* 结束读取计时 */
     read_end = MPI_Wtime();
     read_time = read_end - read_start;
     /* 使用MPI_Reduce收集所有进程的读取时间，取最大值 */
     MPI_Reduce(&read_time, &total_read_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
 
     /* 开始计算计时 */
     compute_start = MPI_Wtime();
 
     /* 计算本地时间平均值 */
     local_avg = (float *)malloc(spatial_size * sizeof(float));
     if (local_avg == NULL) {
         printf("Error: Memory allocation failed for local_avg\n");
         MPI_Abort(MPI_COMM_WORLD, -1);
         return 1;
     }
     
     /* 初始化局部平均值缓冲区 */
     for (i = 0; i < spatial_size; i++) {
         local_avg[i] = 0.0f;
     }
     
     /* 计算局部时间平均值 */
     for (i = 0; i < my_time_count; i++) {
         for (j = 0; j < spatial_size; j++) {
             local_avg[j] += buffer[i * spatial_size + j];
         }
     }
     
     /* 对局部平均值进行归一化 */
     for (j = 0; j < spatial_size; j++) {
         if (my_time_count > 0) {
             local_avg[j] /= my_time_count;
         }
     }
     
     /* 释放原始数据缓冲区，不再需要 */
     free(buffer);
     
     /* 分配全局平均值缓冲区 */
     global_avg = (float *)malloc(spatial_size * sizeof(float));
     if (global_avg == NULL) {
         printf("Error: Memory allocation failed for global_avg\n");
         MPI_Abort(MPI_COMM_WORLD, -1);
         return 1;
     }
     
     // /* 使用MPI归约操作计算全局平均值 - 先求和 */
     // MPI_Reduce(local_avg, global_avg, spatial_size, MPI_FLOAT, MPI_SUM, 0, file_comm);
     
     // /* 进程0对结果进行归一化 */
     // if (proc_in_group == 0) {
     //     for (j = 0; j < spatial_size; j++) {
     //         global_avg[j] /= procs_per_group;
     //     }
     // }
     
     // /* 广播全局平均值给组内所有进程 */
     // MPI_Bcast(global_avg, spatial_size, MPI_FLOAT, 0, file_comm);
 
     /* 使用MPI归约操作计算全局平均值 - 先求和 */
     MPI_Allreduce(local_avg, global_avg, spatial_size, MPI_FLOAT, MPI_SUM, file_comm);
     
     /* 进程0对结果进行归一化 */
 
     for (j = 0; j < spatial_size; j++) {
         global_avg[j] /= procs_per_group;
     }
     
     /* 结束计算计时 */
     compute_end = MPI_Wtime();
     compute_time = compute_end - compute_start;
     /* 使用MPI_Reduce收集所有进程的计算时间，取最大值 */
     MPI_Reduce(&compute_time, &total_compute_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
 
     /* 同步所有进程，确保所有输入文件都已读取和处理 */
     MPI_Barrier(MPI_COMM_WORLD);
     
     if (global_rank == 0) {
         printf("所有输入文件读取和处理完成，开始创建输出文件...\n");
     }
    //  printf("BBBB");
     /* 开始写入计时 */
     write_start_time = MPI_Wtime();
 
     /* 第二阶段：创建输出文件 */
     ret = ncmpi_create(MPI_COMM_WORLD, output_file, NC_64BIT_DATA, info, &ncid_out);
     CHECK_ERR(ret);
     
     /* 为输出文件创建空间维度 (y, x) */
     dimids_out = (int *)malloc(2 * sizeof(int)); /* 只需要2个维度：y和x */
     dim_sizes_out = (MPI_Offset *)malloc(2 * sizeof(MPI_Offset));
     
     dim_sizes_out[0] = dim_sizes_in[1]; /* y 维度大小 */
     dim_sizes_out[1] = dim_sizes_in[2]; /* x 维度大小 */
     
    // printf("CCCC");
     /* 广播 dimids_out 和 dim_sizes_out，确保所有进程值一致 */
     MPI_Bcast(dimids_out, 2, MPI_INT, 0, file_comm);
     MPI_Bcast(dim_sizes_out, 2, MPI_OFFSET, 0, file_comm);
 
     /* 创建y维度 */
     ret = ncmpi_def_dim(ncid_out, dim_names[1], dim_sizes_out[0], &dimids_out[0]);
     CHECK_ERR(ret);
             
     /* 创建x维度 */
     ret = ncmpi_def_dim(ncid_out, dim_names[2], dim_sizes_out[1], &dimids_out[1]);
     CHECK_ERR(ret);
 
     
     /* 创建输出变量 - 使用文件名中的变量类型名作为变量名 */
     /* 为当前文件组定义输出变量 */
     varid_out = (int *)malloc(num_files * sizeof(int));
     /* 使用循环定义每个变量 */
     for (int i = 0; i < num_files; i++) {
         ret = ncmpi_def_var(ncid_out, var_types[i], NC_FLOAT, 2, dimids_out, &varid_out[i]);
         CHECK_ERR(ret);
         /* 添加变量属性，说明这是时间平均值 */
         char attr_text[100];
         sprintf(attr_text, "Time average of %s for %04d-%02d", var_types[i], year, month);
         ret = ncmpi_put_att_text(ncid_out, varid_out[i], "long_name", strlen(attr_text), attr_text);
         CHECK_ERR(ret);
     }
     
     /* 添加全局属性，说明这是时间平均值 */
     char global_attr_text[100];
     sprintf(global_attr_text, "Time average of %s for %04d-%02d", &var_string, year, month);
     ret = ncmpi_put_att_text(ncid_out, NC_GLOBAL, "long_name", strlen(global_attr_text), global_attr_text);
     CHECK_ERR(ret);
 
     /* 结束定义模式 */
     ret = ncmpi_enddef(ncid_out);
     CHECK_ERR(ret);
    //  printf("DDDD");
     
     /* 第三阶段：写入输出文件 */
     if (global_rank == 0) {
         printf("开始写入输出文件...\n");
     }
     
     
     /* 设置写入时按y维度分割的起始位置和计数 */
     MPI_Offset y_size = dim_sizes_out[0];
     MPI_Offset x_size = dim_sizes_out[1];
     MPI_Offset y_chunk = y_size / procs_per_group;
     MPI_Offset y_remainder = y_size % procs_per_group;
     MPI_Offset my_y_count = (proc_in_group < y_remainder) ? y_chunk + 1 : y_chunk;
     MPI_Offset my_y_start = (proc_in_group < y_remainder) ? proc_in_group * (y_chunk + 1) : proc_in_group * y_chunk + y_remainder;
 
     /* 设置写入的起始位置和计数 */
     MPI_Offset write_start[2], write_count[2];
     
     /* 计算该进程处理的数据大小 */
     MPI_Offset proc_data_size = my_y_count * x_size;
     
     /* 创建该进程的数据缓冲区 */
     float *proc_buffer = (float *)malloc(proc_data_size * sizeof(float));
     if (proc_buffer == NULL) {
         printf("Error: Memory allocation failed for proc_buffer\n");
         MPI_Abort(MPI_COMM_WORLD, -1);
         return 1;
     }
     
     /* 复制该进程负责的部分数据 */
     for (i = 0; i < my_y_count; i++) {
         for (j = 0; j < x_size; j++) {
             MPI_Offset global_idx = (my_y_start + i) * x_size + j;
             MPI_Offset local_idx = i * x_size + j;
             proc_buffer[local_idx] = global_avg[global_idx];
         }
     }
    //  printf("EEEE");
     /* 为所有变量写入数据，每个组的进程只实际写入其对应的变量数据 */
     for (int i = 0; i < num_files; i++) {
         if (i == file_group) {
            // printf("c1");
             /* 当前组负责的变量：实际写入数据 */
             write_start[0] = my_y_start;
             write_count[0] = my_y_count;
             write_start[1] = 0;
             write_count[1] = x_size;
            //  printf("c1-1111");
             ret = ncmpi_put_vara_float_all(ncid_out, varid_out[i], write_start, write_count, proc_buffer);
            // ret = ncmpi_put_vara_all(ncid_out, varid_out[i], write_start, write_count, proc_buffer, proc_data_size, MPI_FLOAT);
            //  CHECK_ERR(ret);
         } else {
            // printf("c2");
             /* 其他组负责的变量：count设为0，不实际写入数据 */
             write_start[0] = 0;
             write_count[0] = 0;
             write_start[1] = 0;
             write_count[1] = 0;
            //  printf("c1-2222");
             ret = ncmpi_put_vara_float_all(ncid_out, varid_out[i], write_start, write_count, NULL);
            //  ret = ncmpi_put_vara_all(ncid_out, varid_out[i], write_start, write_count, NC_FLOAT, NULL, 0, NC_FLOAT);
             CHECK_ERR(ret);
             }
     }
     
     /* 关闭输出文件 */
     ret = ncmpi_close(ncid_out);
     CHECK_ERR(ret);
    //  printf("FFFF");
     /* 结束写入计时 */
     write_end_time = MPI_Wtime();
     write_time = write_end_time - write_start_time;
     /* 使用MPI_Reduce收集所有进程的计算时间，取最大值 */
     MPI_Reduce(&write_time, &total_write_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    //  printf("????");
     /* 释放资源 */
     free(dimids_in);
     free(dimids_out);
     free(dim_sizes_in);
     free(dim_sizes_out);
    //  printf("!!!!");
    //  for (i = 0; i < ndims; i++) {
    //      free(dim_names[i]);
    //  }
    free(dim_names);
    //  printf("1");
    //  printf("2");
     free(local_avg);
    //  printf("3");
     free(global_avg);
    //  printf("4");
     free(proc_buffer);
    //  printf("GGGG");
     for (i = 0; i < MAX_FILES; i++) {
         free(input_files[i]);
     }
     free(input_files);
     
     for (i = 0; i < num_var_types; i++) {
         free(var_types[i]);
     }
     free(var_types);
     
     if (proc_in_group == 0) {
         free(varid_out);
     }
    //  printf("HHHH");
     MPI_Info_free(&info);
     MPI_Comm_free(&file_comm);
     
     if (global_rank == 0) {
         printf("成功完成！输出文件: %s\n", output_file);
     }
     /* 结束计时 - 整个程序开始*/
     end_time = MPI_Wtime();
     process_time = end_time - start_time;
     /* 使用MPI_Reduce收集所有进程的读取时间，取最大值 */
     MPI_Reduce(&process_time, &total_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    //  printf("IIII");
     if (global_rank == 0) {
     printf("===== 性能统计 =====\n");
     printf("总读取时间: %.4f 秒\n", total_read_time);
     printf("总计算时间: %.4f 秒\n", total_compute_time);
     printf("总写入时间: %.4f 秒\n", total_write_time);
     printf("总执行时间: %.4f 秒\n", total_time);
     }
     MPI_Finalize();
     return 0;
 }