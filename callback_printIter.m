if mod(iter,100)==0
    fprintf('.');
end
if mod(iter,1000)==0
    fprintf('%0.0f\n',iter);
end