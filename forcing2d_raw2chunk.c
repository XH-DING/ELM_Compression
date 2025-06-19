/*
 * Parallel NetCDF Processor
 * 
 * This program reads a NetCDF file with FLDS/FSDS variable, distributes time steps
 * across processes, and creates a new NetCDF file with all variables and attributes.
 * 
 * Compile with: mpicc -o pnetcdf_processor pnetcdf_processor.c -lpnetcdf
 * Run with: mpiexec -n <num_procs> ./pnetcdf_processor <input_file> <output_file>
 */

 #include <stdio.h>
 #include <stdlib.h>
 #include <string.h>
 #include <pnetcdf.h>
 #include <mpi.h>
 
 #define ERR(e) {if(e) {fprintf(stderr, "Error at %s:%d: %s\n", __FILE__, __LINE__, ncmpi_strerror(e)); MPI_Abort(MPI_COMM_WORLD, 1);}}
 #define MAX_DIMS 10
 #define MAX_VARS 100
 #define MAX_ATTR_NAME 256
 #define MAX_ATTR_VAL 10240
 #define MAX_VAR_NAME 256
 
 static MPI_Datatype
nc2mpitype(nc_type xtype)
{
    switch(xtype){
        case NC_CHAR :   return MPI_CHAR;
        case NC_BYTE :   return MPI_SIGNED_CHAR;
        case NC_SHORT :  return MPI_SHORT;
        case NC_INT :    return MPI_INT;
        case NC_FLOAT :  return MPI_FLOAT;
        case NC_DOUBLE : return MPI_DOUBLE;
        case NC_UBYTE :  return MPI_UNSIGNED_CHAR;
        case NC_USHORT : return MPI_UNSIGNED_SHORT;
        case NC_UINT :   return MPI_UNSIGNED;
        case NC_INT64 :  return MPI_LONG_LONG_INT;
        case NC_UINT64 : return MPI_UNSIGNED_LONG_LONG;
        default:         return MPI_DATATYPE_NULL;
    }
}

 int main(int argc, char *argv[]) {
     int rank, nprocs, ret;
     MPI_Init(&argc, &argv);
     MPI_Comm_rank(MPI_COMM_WORLD, &rank);
     MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
     
    //  volatile int attached = 0;
    //  while (!attached) {
    //     sleep(1); // Wait for all processes to attach
    //     if (getenv("DEBUG_READY")) attached = 1;
    //  }
     if (argc != 3) {
         if (rank == 0) {
             printf("Usage: %s <input_file> <output_file>\n", argv[0]);
         }
         MPI_Finalize();
         return 1;
     }
 
     char *input_file = argv[1];
     char *output_file = argv[2];
     
     int ncid_in, ncid_out;
     int ndims, nvars, natts, unlimdimid;
     char main_var_name[MAX_VAR_NAME];
     
     // Get the main variable name from the input file path
     char *base_name = strrchr(input_file, '/');
     if (base_name == NULL) base_name = input_file;
     else base_name++;
     
    // 处理文件名，提取变量名
     // 拷贝 filename（strtok 需要修改字符串）
     char name_copy[MAX_VAR_NAME];
     strncpy(name_copy, base_name, MAX_VAR_NAME);
     name_copy[MAX_VAR_NAME - 1] = '\0';
 
     // 分割 filename
     char *token;
     char *tokens[10];  // 最多支持10段
     int token_count = 0;
 
     token = strtok(name_copy, ".");
     while (token != NULL && token_count < 10) {
         tokens[token_count++] = token;
         token = strtok(NULL, ".");
     }
 
     if (token_count >= 4) {
         strncpy(main_var_name, tokens[3], MAX_VAR_NAME);
         main_var_name[MAX_VAR_NAME - 1] = '\0';
         if (rank == 0) printf("Detected main variable: %s\n", main_var_name);
     } else {
         if (rank == 0) {
             printf("Failed to extract variable name from file: %s\n", base_name);
             printf("Expected format: clmforc.Daymet4.1km.VARNAME.YYYY-MM.nc\n");
         }
         MPI_Finalize();
         return 1;
     }

     if (rank == 0) {
         printf("Processing file: %s\n", input_file);
         printf("Main variable: %s\n", main_var_name);
         printf("Output file: %s\n", output_file);
     }
     
     printf("1111\n");
     // Open input file
     ret = ncmpi_open(MPI_COMM_WORLD, input_file, NC_NOWRITE, MPI_INFO_NULL, &ncid_in);
     ERR(ret);
    //  printf("****\n");
     // Get file information
     ret = ncmpi_inq(ncid_in, &ndims, &nvars, &natts, &unlimdimid);
     ERR(ret);
     
     if (rank == 0) {
         printf("Input file has %d dimensions, %d variables, %d global attributes\n", 
                ndims, nvars, natts);
     }
     
     // Get dimension information
     int dim_ids[MAX_DIMS];
     char dim_names[MAX_DIMS][MAX_VAR_NAME];
     MPI_Offset dim_lens[MAX_DIMS];
     
     for (int i = 0; i < ndims; i++) {
         ret = ncmpi_inq_dim(ncid_in, i, dim_names[i], &dim_lens[i]);
         ERR(ret);
         dim_ids[i] = i;
         if (rank == 0) {
             printf("Dimension %d: %s, length = %lld\n", i, dim_names[i], dim_lens[i]);
         }
     }
     
     // Find time dimension
     int time_dim_id = -1;
     for (int i = 0; i < ndims; i++) {
         if (strcmp(dim_names[i], "time") == 0) {
             time_dim_id = i;
             break;
         }
     }
     
     if (time_dim_id == -1) {
         if (rank == 0) {
             printf("Error: Cannot find 'time' dimension in the input file.\n");
         }
         ncmpi_close(ncid_in);
         MPI_Finalize();
         return 1;
     }

     MPI_Info info;
     MPI_Info_create(&info);
     MPI_Info_set(info, "nc_chunk_default_filter", "sz");
     MPI_Info_set(info, "nc_chunking", "enable");
     // Create output file
     ret = ncmpi_create(MPI_COMM_WORLD, output_file, NC_CLOBBER, info, &ncid_out);
     ERR(ret);
    //  printf("2222\n");
     // Define dimensions in output file
     int out_dim_ids[MAX_DIMS];
     for (int i = 0; i < ndims; i++) {
         if (i == unlimdimid) {
            ret = ncmpi_def_dim(ncid_out, dim_names[i], NC_UNLIMITED, &out_dim_ids[i]);
            // ret = ncmpi_def_dim(ncid_out, dim_names[i], dim_lens[i], &out_dim_ids[i]);
         } else {
             ret = ncmpi_def_dim(ncid_out, dim_names[i], dim_lens[i], &out_dim_ids[i]);
         }
         ERR(ret);
     }
     
     // Copy global attributes
     for (int i = 0; i < natts; i++) {
         char att_name[MAX_ATTR_NAME];
         nc_type att_type;
         MPI_Offset att_len;
         
         ret = ncmpi_inq_attname(ncid_in, NC_GLOBAL, i, att_name);
         ERR(ret);
         
         ret = ncmpi_inq_att(ncid_in, NC_GLOBAL, att_name, &att_type, &att_len);
         ERR(ret);
         
         void *att_val = malloc(att_len * sizeof(char) * MAX_ATTR_VAL);
         ret = ncmpi_get_att(ncid_in, NC_GLOBAL, att_name, att_val);
         ERR(ret);
         
         ret = ncmpi_put_att(ncid_out, NC_GLOBAL, att_name, att_type, att_len, att_val);
         ERR(ret);
         
         free(att_val);
     }

     // Define variables in output file
     int var_ids[MAX_VARS];
     int out_var_ids[MAX_VARS];
     char var_names[MAX_VARS][MAX_VAR_NAME];
     int main_var_id = -1;
     int out_main_var_id = -1;
     
     for (int i = 0; i < nvars; i++) {
         int var_ndims;
         int var_dimids[MAX_DIMS];
         int var_natts;
         nc_type var_type;
         
         ret = ncmpi_inq_var(ncid_in, i, var_names[i], &var_type, &var_ndims, var_dimids, &var_natts);
         ERR(ret);
         
         var_ids[i] = i;
         
         // Map input dimension IDs to output dimension IDs
         int out_var_dimids[MAX_DIMS];
         for (int j = 0; j < var_ndims; j++) {
             out_var_dimids[j] = out_dim_ids[var_dimids[j]];
         }
         
         ret = ncmpi_def_var(ncid_out, var_names[i], var_type, var_ndims, out_var_dimids, &out_var_ids[i]);
         ERR(ret);
         
         // If this is the main variable, save its ID
         if (strcmp(var_names[i], main_var_name) == 0) {
             main_var_id = i;
             out_main_var_id = out_var_ids[i];
             int chunk_dim[3] = {1, dim_lens[1], dim_lens[2]};
             ret = ncmpi_var_set_chunk(ncid_out, out_main_var_id, chunk_dim);
             ERR(ret);
            //  ret = ncmpi_var_set_filter(ncid_out, out_main_var_id, NC_FILTER_SZ);
            //  ERR(ret);
             if (rank == 0) {
                 printf("Found main variable %s with ID %d\n", main_var_name, main_var_id);
             }
         }
         else {
            continue; // Skip other variables for now
            // int chunk_dim[var_ndims];
            // // 获取每个维度长度
            // for (int d = 0; d < var_ndims; d++) {
            //     chunk_dim[d] = dim_lens[var_dimids[d]];
            // }
            // // 设置 chunk 大小
            // // ret = ncmpi_var_set_chunk(ncid_out, out_var_ids[i], chunk_dim);
            // // ret = ncmpi_var_set_filter(ncid_out, out_var_ids[i], NC_FILTER_NONE);
            //  if (rank == 0) {
            //      printf("Found variable %s with ID %d\n", var_names[i], i);
            //  }
         }
         
         // Copy variable attributes
         for (int j = 0; j < var_natts; j++) {
             char att_name[MAX_ATTR_NAME];
             nc_type att_type;
             MPI_Offset att_len;
             
             ret = ncmpi_inq_attname(ncid_in, i, j, att_name);
             ERR(ret);
             
             ret = ncmpi_inq_att(ncid_in, i, att_name, &att_type, &att_len);
             ERR(ret);
             
             void *att_val = malloc(att_len * sizeof(char) * MAX_ATTR_VAL);
             ret = ncmpi_get_att(ncid_in, i, att_name, att_val);
             ERR(ret);
             
             ret = ncmpi_put_att(ncid_out, out_var_ids[i], att_name, att_type, att_len, att_val);
             ERR(ret);
             
             free(att_val);
         }
     }
     
     if (main_var_id == -1) {
         if (rank == 0) {
             printf("Error: Cannot find the main variable '%s' in the input file.\n", main_var_name);
         }
         ncmpi_close(ncid_in);
         ncmpi_close(ncid_out);
         MPI_Finalize();
         return 1;
     }
     
     // End define mode for output file
     ret = ncmpi_enddef(ncid_out);
     ERR(ret);
     printf("3333\n");
     // Process the main variable by distributing time steps
     int main_var_ndims;
     int main_var_dimids[MAX_DIMS];
     nc_type main_var_type;
     
     ret = ncmpi_inq_var(ncid_in, main_var_id, NULL, &main_var_type, &main_var_ndims, main_var_dimids, NULL);
     ERR(ret);
    //  printf("type %d\n", main_var_type);

     // Find the time dimension in the main variable
     int time_dim_index = -1;
     for (int i = 0; i < main_var_ndims; i++) {
         if (main_var_dimids[i] == time_dim_id) {
             time_dim_index = i;
             break;
         }
     }
     
     if (time_dim_index == -1) {
         if (rank == 0) {
             printf("Error: Main variable does not have the time dimension.\n");
         }
         ncmpi_close(ncid_in);
         ncmpi_close(ncid_out);
         MPI_Finalize();
         return 1;
     }
     
     // Calculate time steps distribution per process
     MPI_Offset time_len = dim_lens[time_dim_id];
     MPI_Offset time_per_proc = time_len / nprocs;
     MPI_Offset time_remainder = time_len % nprocs;
     
     MPI_Offset start_time = rank * time_per_proc + (rank < time_remainder ? rank : time_remainder);
     MPI_Offset count_time = time_per_proc + (rank < time_remainder ? 1 : 0);
     
     if (rank == 0) {
         printf("Total time steps: %lld\n", time_len);
         printf("Distributing across %d processes\n", nprocs);
     }
     
     printf("Process %d: processing time steps %lld to %lld (count: %lld)\n", 
            rank, start_time, start_time + count_time - 1, count_time);
     
     // Prepare start and count arrays for reading/writing
     MPI_Offset start[MAX_DIMS], count[MAX_DIMS];
     for (int i = 0; i < main_var_ndims; i++) {
         if (i == time_dim_index) {
             start[i] = start_time;
             count[i] = count_time;
         } else {
             start[i] = 0;
             count[i] = dim_lens[main_var_dimids[i]];
         }
        //  printf("start[%d] = %lld, count[%d] = %lld\n", i, start[i], i, count[i]);
     }
     
     // Calculate buffer size and allocate memory
     MPI_Offset buffer_size = 1;
     for (int i = 0; i < main_var_ndims; i++) {
         buffer_size *= count[i];
     }
     
     void *buffer = NULL;
     switch (main_var_type) {
         case NC_BYTE:
         case NC_CHAR:
             buffer = malloc(buffer_size * sizeof(char));
             break;
         case NC_SHORT:
             buffer = malloc(buffer_size * sizeof(short));
             break;
         case NC_INT:
             buffer = malloc(buffer_size * sizeof(int));
             break;
         case NC_FLOAT:
             buffer = malloc(buffer_size * sizeof(float));
             break;
         case NC_DOUBLE:
             buffer = malloc(buffer_size * sizeof(double));
             break;
         default:
             if (rank == 0) {
                 printf("Error: Unsupported variable type %d\n", main_var_type);
             }
             ncmpi_close(ncid_in);
             ncmpi_close(ncid_out);
             MPI_Finalize();
             return 1;
     }

     if (buffer == NULL) {
         printf("Error: Failed to allocate buffer of size %lld bytes on process %d\n", 
                buffer_size * sizeof(double), rank);
         ncmpi_close(ncid_in);
         ncmpi_close(ncid_out);
         MPI_Finalize();
         return 1;
     }

     // Read the main variable
     ret = ncmpi_get_vara_all(ncid_in, main_var_id, start, count, buffer, buffer_size, nc2mpitype(main_var_type));
     ERR(ret);

    //  // 打印变量 ID 和类型
    // printf("[Rank %d] Calling ncmpi_put_vara_all:\n", rank);
    // printf("  -> ncid_out        = %d\n", ncid_out);
    // printf("  -> out_main_var_id = %d\n", out_main_var_id);
    // printf("  -> main_var_type   = %d (MPI type = %d)\n", main_var_type, nc2mpitype(main_var_type));

    // // 打印 start 和 count 数组内容
    // printf("  -> start = [");
    // for (int i = 0; i < main_var_ndims; i++) {
    //     printf("%lld", start[i]);
    //     if (i != main_var_ndims - 1) printf(", ");
    // }
    // printf("]\n");

    // printf("  -> count = [");
    // for (int i = 0; i < main_var_ndims; i++) {
    //     printf("%lld", count[i]);
    //     if (i != main_var_ndims - 1) printf(", ");
    // }
    // printf("]\n");

    // // 打印 buffer 地址和大小
    // printf("  -> buffer = %p\n", buffer);
    // printf("  -> buffer_size = %lld (elements)\n", buffer_size);

    // float new_buffer[4][8075][7814];
    // int i,j,k;
    // for (i=0; i<8075; i++)
    //     for (j=0; j<7814; j++)
    //         for (k=0; k<4; k++)
    //             new_buffer[k][i][j] = 0.0;
    double write_start_time, write_end_time, write_time, total_write_time;
    write_start_time = MPI_Wtime();
    printf("vvvvv\n");
     ret = ncmpi_put_vara_all(ncid_out, out_main_var_id, start, count, buffer, buffer_size, nc2mpitype(main_var_type));
     ERR(ret);
     printf("wwwww\n");
    //  ret = ncmpi_put_vara_float_all(ncid_out, out_main_var_id, start, count, buffer);
    //  ERR(ret);
    write_end_time = MPI_Wtime();
    write_time = write_end_time - write_start_time;
    /* 使用MPI_Reduce收集所有进程的计算时间，取最大值 */
    MPI_Reduce(&write_time, &total_write_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
     if (rank == 0) {
        printf("总写入时间: %.4f 秒\n", total_write_time);
     }
     printf("Rank: %d, Write Time: %.4f\n", rank, write_time);
     free(buffer);
     
    //  // Copy other variables
    //  for (int i = 0; i < nvars; i++) {
    //      if (i == main_var_id) {
    //          continue; // Skip the main variable (already processed)
    //      }
         
    //      // Get variable info
    //      int var_ndims;
    //      int var_dimids[MAX_DIMS];
    //      nc_type var_type;
         
    //      ret = ncmpi_inq_var(ncid_in, i, NULL, &var_type, &var_ndims, var_dimids, NULL);
    //      ERR(ret);
         
    //      // Check if variable has time dimension
    //      int has_time_dim = 0;
    //      int time_index = -1;
    //      for (int j = 0; j < var_ndims; j++) {
    //          if (var_dimids[j] == time_dim_id) {
    //              has_time_dim = 1;
    //              time_index = j;
    //              break;
    //          }
    //      }
         
    //      // Prepare start and count arrays
    //      MPI_Offset var_start[MAX_DIMS], var_count[MAX_DIMS];
    //      for (int j = 0; j < var_ndims; j++) {
    //          if (has_time_dim && j == time_index) {
    //              var_start[j] = start_time;
    //              var_count[j] = count_time;
    //          } else {
    //              var_start[j] = 0;
    //              var_count[j] = dim_lens[var_dimids[j]];
    //          }
    //      }
         
    //      // Calculate buffer size
    //      MPI_Offset var_buffer_size = 1;
    //      for (int j = 0; j < var_ndims; j++) {
    //          var_buffer_size *= var_count[j];
    //      }

    //      // Allocate buffer
    //      void *var_buffer = NULL;
    //      switch (var_type) {
    //          case NC_BYTE:
    //          case NC_CHAR:
    //              var_buffer = malloc(var_buffer_size * sizeof(char));
    //              break;
    //          case NC_SHORT:
    //              var_buffer = malloc(var_buffer_size * sizeof(short));
    //              break;
    //          case NC_INT:
    //              var_buffer = malloc(var_buffer_size * sizeof(int));
    //              break;
    //          case NC_FLOAT:
    //              var_buffer = malloc(var_buffer_size * sizeof(float));
    //              break;
    //          case NC_DOUBLE:
    //              var_buffer = malloc(var_buffer_size * sizeof(double));
    //              break;
    //          default:
    //              if (rank == 0) {
    //                  printf("Error: Unsupported variable type %d for variable %s\n", var_type, var_names[i]);
    //              }
    //              continue;
    //      }
         
    //      if (var_buffer == NULL) {
    //          printf("Error: Failed to allocate buffer of size %lld bytes for variable %s on process %d\n", 
    //                 var_buffer_size * sizeof(double), var_names[i], rank);
    //          continue;
    //      }
         
    //      // If variable has time dimension, each process reads/writes its portion
    //      if (has_time_dim) {
    //          ret = ncmpi_get_vara_all(ncid_in, i, var_start, var_count, var_buffer, var_buffer_size, nc2mpitype(var_type));
    //          ERR(ret);
             
    //         //  ret = ncmpi_put_vara_all(ncid_out, out_var_ids[i], var_start, var_count, var_buffer, var_buffer_size, nc2mpitype(var_type));
    //         //  ERR(ret);
    //      } 
    //      // If variable doesn't have time dimension, only process 0 reads/writes it
    //      else {
            
    //          ret = ncmpi_get_var_all(ncid_in, i, var_buffer, var_buffer_size, nc2mpitype(var_type));
    //          ERR(ret);   

    //         //  ret = ncmpi_put_var_all(ncid_out, out_var_ids[i], var_buffer, var_buffer_size, nc2mpitype(var_type));
    //         //  ERR(ret);   
    //         }
    //     //  else if (rank == 0) {
    //     //      ret = ncmpi_get_vara(ncid_in, i, var_start, var_count, var_buffer, var_buffer_size, nc2mpitype(var_type));
    //     //      ERR(ret);
             
    //         //  ret = ncmpi_put_vara(ncid_out, out_var_ids[i], var_start, var_count, var_buffer, var_buffer_size, nc2mpitype(var_type));
    //         //  ERR(ret);
    //     //  }
         
    //      free(var_buffer);
    //  }
     
     // Close files
     ret = ncmpi_close(ncid_in);
     ERR(ret);
     
     ret = ncmpi_close(ncid_out);
     ERR(ret);
    //  MPI_Info_free(&info);
     
     if (rank == 0) {
         printf("Processing completed successfully.\n");
     }
     
     MPI_Finalize();
     return 0;
 }
