testNwaymLike = function(g,y,Xm)
{
	P = dim(Xm)[2]
	N = dim(Xm)[1]
	
	y = matrix(y,N)
	Xc = t(t(Xm) - colMeans(Xm))
	
	m0 = -0.5*(N - 1) * log(var(y)*(N-1))
	
	.Call("RjeffmlikeNWayAov",
		t(Xc)%*%Xc,
		t(Xc)%*%y,
		var(y)*(N-1),
		N, P, g, 
		package="BayesFactor") - m0
}

