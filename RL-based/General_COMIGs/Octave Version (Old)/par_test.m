function retval = par_test (num, fun)
  tic;
  result = pararrayfun(nproc, fun, num);
  retval = toc;
endfunction
