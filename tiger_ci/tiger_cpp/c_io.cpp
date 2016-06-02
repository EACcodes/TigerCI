/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */
#include "c_io.h"
#include "StringUtils.h"

// Global table of file names, ugly but only used in this file
unordered_map<string, FILE*> c_io_files;

std::string file_name_clean(const char* filename){
  std::string fixed(filename);
  return trim(fixed);
}

void c_close_files(){
  for(auto item = c_io_files.begin(); item != c_io_files.end(); ++item){
    FILE* fp = item->second;
    int ret_code = fclose(fp);
    if (ret_code != 0){
      cout << "Error code = " << ret_code << " while trying to close file " << item->first << endl;
    }
  }
}
    
void c_setup_vector_io(const char* filename){
  std::string fname = file_name_clean(filename);
  FILE* fp = fopen(fname.c_str(), "wb+");
  if (!fp) {
    cout << "Hit an error trying to open file " << filename << endl;
    cout << "error code = " << ferror(fp) << endl;
  }
  string name {fname};
  c_io_files[fname] = fp;
}

FILE* get_file_pointer(const char* fname){  
  std::string filename = file_name_clean(fname);
  auto search = c_io_files.find(filename);
  if (search == c_io_files.end()){
    cout << "Fatal error ... tried to access file " << filename << " with out running the setup_vector_io routine first!" << endl;
    exit(2);
  }
  return search->second;
}

void c_write_vector(const char* filename, int64_t& position, double* vector, int64_t& length){
  FILE* fp = get_file_pointer(filename);
  fseek(fp, (position-1)*length*sizeof(double), SEEK_SET);
  size_t count = fwrite(vector, sizeof(double), length , fp);
  error_processing(fp, filename, "writing", count, length, position);  
}

void c_read_vector(const char* filename, int64_t& position, double* vector, int64_t& length){
   FILE* fp = get_file_pointer(filename);
  fseek(fp, (position-1)*length*sizeof(double), SEEK_SET);
  size_t count = fread(vector, sizeof(double), length, fp);
  error_processing(fp, filename, "reading", count, length, position);
}

void c_write_ints(const char* filename, int64_t& position, int64_t* vector, int64_t& length){
  FILE* fp = get_file_pointer(filename);
  fseek(fp, (position-1)*length*sizeof(int64_t), SEEK_SET);
  size_t count = fwrite(vector, sizeof(int64_t), length , fp);
  error_processing(fp, filename, "writing", count, length, position);  
}

void c_read_ints(const char* filename, int64_t& position, int64_t* vector, int64_t& length){
   FILE* fp = get_file_pointer(filename);
  fseek(fp, (position-1)*length*sizeof(int64_t), SEEK_SET);
  size_t count = fread(vector, sizeof(int64_t), length, fp);
  error_processing(fp, filename, "reading", count, length, position);
}


void error_processing(FILE* fp, const char* filename, const char* rw, size_t count, size_t length, size_t position){
  if (!(count==length)){
    cout << endl << endl;
    cout << "   Encountered an error while " << rw << " from the file " << filename << endl;
    cout << "   "<<rw << " vector at position " << position << endl;
    cout << "   found " << count << " numbers but " << length << " were requested" << endl;
    if (feof(fp)){
      cout << "   error occured because the end of the file was reached " << endl;
    }
    else if (ferror(fp)){
      cout << "  some error occured: ferror = " << ferror(fp) << endl;
    }
    else {
      cout << "  something really weird happened (both ferror and feof are zero?)" << endl;
    }
    cout  << endl << endl ; 
    exit(2);
  }
}
