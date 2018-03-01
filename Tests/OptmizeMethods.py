import scipy.optimize as opt


def func(x):
    y = (x+1)**2
    return y


MyRan = slice(-5, 5, 2)
InitGuess = opt.brute(func, (MyRan,))
result = opt.minimize(func, InitGuess, method="BFGS")

MinPoint = result.get("x")
MinValue = result.get("fun")
print("Valor de x: " + str(MinPoint))
print("Valor minimo: " + str(MinValue))

