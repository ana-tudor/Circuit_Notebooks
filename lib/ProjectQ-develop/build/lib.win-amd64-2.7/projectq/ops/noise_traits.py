#-----
import random, math

scale = 3E-6

Rx_pdf       = random.gauss
Rx_pdf_args  = 0., scale*math.pi

Rz_pdf       = random.gauss
Rz_pdf_args  = 0., scale*math.pi

def normal2(sigma1, sigma2):
     return random.gauss(0., sigma1), random.gauss(0., sigma2)

CNOT_pdf       = normal2
CNOT_pdf_args  = scale*math.pi, scale*math.pi
#CNOT_epsilon   = 0.
#-----
