function y = array_response(a1,a2,N)
for m= 0:sqrt(N)-1
    for n= 0:sqrt(N)-1
        y(m*(sqrt(N))+n+1) = exp( 1i* pi* ( m*sin(a1)*sin(a2) + n*cos(a2) ) );
    end
end
y = y.'/sqrt(N);
end