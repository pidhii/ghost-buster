let g:ale_cpp_gcc_options = system("root-config --cflags")[0:-2] . " -I include -std=gnu++17 " 
