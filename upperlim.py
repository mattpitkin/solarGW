# Getting an upper limit on solarGWs
def upperlim(p,duration):
	Xspacing = 2.44140625E-4
	N = duration/Xspacing
	h = sqrt(p/N)
	print h
