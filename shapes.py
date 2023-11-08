import numpy as np
import numba as nb


@nb.njit(inline="always")
def integrand(x, y, params):
    # TODO: Replace placeholder integrand
    return params[0] + params[1]*x*x + params[2]*y*y + params[3]*np.sin(x) + params[4]*np.cos(y)


@nb.njit(inline="always")
def dunavant7(newF, x1, y1, f1, x2, y2, f2, x3, y3, f3, params):
    total = (f1+f2+f3)/40. + newF/15.
    total += integrand((x1+x2)/2., (y1+y2)/2, params) / 15.
    total += integrand((x1+x3)/2., (y1+y3)/2, params) / 15.
    total += integrand((x1+x2+x3)/3., (y1+y2+y3)/3, params) * (9./40.)
    jacobian = (x2-x1)*(y3-y1) - (y2-y1)*(x3-x1)
    return total * abs(jacobian)


@nb.njit
def integrateTriangle(stack, params, tol=1e-3, initEval=True):
    '''Stack is a 16 array, with entry 0 containing the right triangle to integrate,
    defined by 3 x-y pairs and their function values starting with the right angled one.'''
    # xyf: 0,1,2     3,4,5    6,7,8
    
    if initEval:
        stack[0,2] = integrand(stack[0,0], stack[0,1], params)
        stack[0,5] = integrand(stack[0,3], stack[0,4], params)
        stack[0,8] = integrand(stack[0,6], stack[0,7], params)

    i = 0
    total = 0.
    while i >= 0:
        x0, y0, f0, x1, y1, f1, x2, y2, f2 = stack[i]
        newX = (x1+x2)/2.
        newY = (y1+y2)/2.
        newF = integrand(newX, newY, params)
        if abs(newF - (f2+f1)/2.) < tol or i >= 12:
            total += dunavant7(newF, x0, y0, f0, x1, y1, f1, x2, y2, f2, params)
            i -= 1
        else:
            # Append two triangles subdivided from the first
            stack[i+1] = newX, newY, newF, x0, y0, f0, x1, y1, f1
            stack[i] = newX, newY, newF, x0, y0, f0, x2, y2, f2
            i += 1
            
    return total
