function a = test_batch(sz,fname)
a = rand(sz,sz)
    
fid = fopen(fname,'w');
fprintf(fid,'%g', a);
end