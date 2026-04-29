# importing libraries 
import matplotlib.pyplot as plt 
import numpy as np 
import sympy as sym 
  
# creating two array for plotting 
X = ['Oct 21', 'Nov 21', 'Dec 21', 'Jan 22', 'Feb 22', 'Mar 22', 'Apr 22', 
     'May 22', 'Jun 22', 'Jul 22', 'Aug 22', 'Sept 22', 'Oct 22', 'Nov 22',
     'Dec 22', 'Jan 23', 'Feb 23', 'Mar 23', 'Apr 23', 'May 23', 'Jun 23',
     'Jul 23', 'Aug 23', 'Sept 23', 'Oct 23']
Y = [1764096, 2189474, 3463682, 4465010, 2218769, 2626130, 1505382, 1979853, 
     1636089, 2636005, 2055372, 2439880, 4278833, 2526506, 3611316, 743170, 
     1922541, 1705668, 1326187, 1029268, 1433912, 880781, 1100544, 1191511, 1733155] 
  
# creating scatter plot with both negative  
# and positive axes 
plt.scatter(X, Y) 
plt.xticks(rotation=90)
  
# adding vertical line in data co-ordinates 
plt.axvline('Oct 22', c='black', ls='--') 
  
# adding horizontal line in data co-ordinates 
plt.axhline(0, c='black', ls='-') 
  
# visualizing the plot using plt.show() function 
plt.annotate('Month of Scandal', xy=(12,4278833), xytext=(15,4278833),
             fontsize=9, arrowprops=dict(arrowstyle='->'))
plt.savefig("RawData.svg")
plt.show()

#%% DD
xs = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24]

x=sym.symbols('x')
def myDD(X,y):
    n=len(X)
    DD=np.zeros([n,n+1]) 
    DD[:,0]=X
    DD[[range(0,n)],1]=y 
    for j in range(0,n-1): 
        DD[j,2]=(DD[j+1,1]-DD[j,1])/(DD[j+1,0]-DD[j,0]) 
    for j in range(3,n+1): 
        for i in range(0,(n+1)-j): 
            DD[i,j]=(DD[i+1,j-1]-DD[i,j-1])/(DD[i+(j-1),0]-DD[i,0]) 
    for i in range(0,n):
        for j in range(0,n+1):
            if i+j>n:
                DD[i,j]=np.nan 
    NIP=0
    for i in range(1,n+1):
        prod=DD[0,i]
        for j in range(0,i-1):
            prod=prod*(x-X[j])
        NIP=NIP+prod
    p=sym.simplify(NIP)
    return(p)

poly=myDD(xs,Y)
print(poly) 

#%% split DD

X1 = np.array([0,1,2,3,4,5,6,7,8,9,10,11,12])
Y1 = [1764096, 2189474, 3463682, 4465010, 2218769, 2626130, 1505382, 1979853, 
     1636089, 2636005, 2055372, 2439880, 4278833]

X2 = np.array([12,13,14,15,16,17,18,19,20,21,22,23,24])
Y2 = [4278833, 2526506, 3611316, 743170, 
     1922541, 1705668, 1326187, 1029268, 1433912, 880781, 1100544, 1191511, 1733155]

poly1=myDD(X1,Y1)
print(poly1)
poly2=myDD(X2,Y2)
print(poly2) 

#%% Cubic Spline
X1 = np.array([0,1,2,3,4,5,6,7,8,9,10,11,12])
Y1 = [1764096, 2189474, 3463682, 4465010, 2218769, 2626130, 1505382, 1979853, 
     1636089, 2636005, 2055372, 2439880, 4278833]

X2 = np.array([12,13,14,15,16,17,18,19,20,21,22,23,24])
Y2 = [4278833, 2526506, 3611316, 743170, 
     1922541, 1705668, 1326187, 1029268, 1433912, 880781, 1100544, 1191511, 1733155]

def myCubicSpline(xs,ys):
    n=len(xs)
    x,a,b,c,d=sym.symbols('x a b c d')
    a=sym.symbols('a:%d'%(n-1))
    b=sym.symbols('b:%d'%(n-1))
    c=sym.symbols('c:%d'%(n-1))
    d=sym.symbols('d:%d'%(n-1))

    sn=np.array([])
    for i in range(0,n-1):
        sn=np.append(sn,a[i]+b[i]*x+c[i]*x**2+d[i]*x**3)

    equations=np.array([])
    for i in range(0,n-1):
        #evaluating sj at t_j and t_j+1
        sjtj=(sn[i]-ys[i]).subs(x,xs[i])
        sjtj1=(sn[i]-ys[i+1]).subs(x,xs[i+1])
        if i<(n-2):
            #smoothness
            smooth=sym.diff(sn[i],x,1).subs(x,xs[i+1])-sym.diff(sn[i+1],x,1).subs(x,xs[i+1])
            #double smoothness
            doublesmooth=sym.diff(sn[i],x,2).subs(x,xs[i+1])-sym.diff(sn[i+1],x,2).subs(x,xs[i+1])
        equations=np.append(equations,[sjtj,sjtj1,smooth,doublesmooth])
    #natural conditions
    equations=np.append(equations,sym.diff(sn[0],x,2).subs(x,xs[0]))
    equations=np.append(equations,sym.diff(sn[-1],x,2).subs(x,xs[-1]))

    coeffs=sym.solve(equations,dict=True)
    sol_array=np.array([])
    for i in range(0,n-1):
        sol_array=np.append(sol_array,coeffs[0][a[i]]+coeffs[0][b[i]]*x+coeffs[0][c[i]]*x**2+coeffs[0][d[i]]*x**3)
    return(sol_array)

spline_funcs=myCubicSpline(xs,Y)
spline_funcs1=myCubicSpline(X1,Y1)
spline_funcs2=myCubicSpline(X2,Y2)
print(sym.simplify(spline_funcs1))
print(sym.simplify(spline_funcs2))
#%% plotting cubic spline

x=sym.symbols('x')
func=0*x
func=sym.Piecewise((func,False),(spline_funcs[0],x<=xs[1]))

for i in range(1,len(xs)-2):
    func=sym.Piecewise((func,x<=xs[i]),(spline_funcs[i],x<=xs[i]+1))

func=sym.Piecewise((func,x<=xs[-1]-1),(spline_funcs[-1],x>xs[-1]-1))

x_new=np.linspace(0,24,500)
y_new=[func.subs(x,Y) for Y in x_new]

plt.plot(x_new,y_new)
plt.xticks(rotation=90)
plt.axhline(0, c='black', ls='-') 
plt.scatter(X, Y, c='darkgreen') 
plt.axvline('Oct 22', c='black', ls='--') 
plt.title('Cubic Spline Interpolation')
plt.xlabel('Month & Year')
plt.ylabel('Average Views Per Month')
plt.savefig("CubicSpline.svg")
plt.show()
#%% plotting filled cubic spline

x=sym.symbols('x')
func=0*x
func=sym.Piecewise((func,False),(spline_funcs[0],x<=xs[1]))

for i in range(1,len(xs)-2):
    func=sym.Piecewise((func,x<=xs[i]),(spline_funcs[i],x<=xs[i]+1))

func=sym.Piecewise((func,x<=xs[-1]-1),(spline_funcs[-1],x>xs[-1]-1))

x_new=np.linspace(0,24,500)
y_new=[func.subs(x,Y) for Y in x_new]

plt.plot(x_new,y_new)
plt.fill_between(x_new.astype(np.float64),np.asarray(y_new).astype(np.float64),
                 where=(0<x_new)&(x_new<12),color='forestgreen',alpha=0.2)
plt.fill_between(x_new.astype(np.float64),np.asarray(y_new).astype(np.float64),
                 where=(12<x_new)&(x_new<24),color='deeppink',alpha=0.1)
plt.annotate('Area=30,081,881.28',xy=(2,700000),xytext=(2,700000),fontsize=10)
plt.annotate('Area=20,179,589.12',xy=(15,300000),xytext=(15,300000),fontsize=10)
plt.xticks(rotation=90)
plt.axhline(0, c='black', ls='-') 
plt.scatter(X, Y, c='darkgreen') 
plt.axvline('Oct 22', c='black', ls='--') 
plt.title('Cubic Spline Interpolation')
plt.xlabel('Month & Year')
plt.ylabel('Average Views Per Month')
plt.savefig("FilledCubicSpline.svg")
plt.show()
#%% plotting cubic spline pre

x1 = ['Oct 21', 'Nov 21', 'Dec 21', 'Jan 22', 'Feb 22', 'Mar 22', 'Apr 22', 
     'May 22', 'Jun 22', 'Jul 22', 'Aug 22', 'Sept 22', 'Oct 22']

x=sym.symbols('x')
func=0*x
func=sym.Piecewise((func,False),(spline_funcs1[0],x<=X1[1]))

for i in range(1,len(X1)-2):
    func=sym.Piecewise((func,x<=X1[i]),(spline_funcs1[i],x<=X1[i]+1))

func=sym.Piecewise((func,x<=X1[-1]-1),(spline_funcs1[-1],x>X1[-1]-1))

x_new=np.linspace(0,12,500)
y_new=[func.subs(x,Y1) for Y1 in x_new]

plt.plot(x_new,y_new)
plt.xticks(rotation=90)
plt.axhline(0, c='black', ls='-') 
plt.scatter(x1, Y1, c='darkgreen') 
plt.title('Cubic Spline Interpolation Pre Scandal')
plt.xlabel('Month & Year')
plt.ylabel('Average Views Per Month')
plt.savefig("PreCubicSpline.svg")
plt.show()
#%% plotting cubic spline post

x2 = ['Oct 22', 'Nov 22', 'Dec 22', 'Jan 23', 'Feb 23', 'Mar 23', 'Apr 23', 
     'May 23', 'Jun 23', 'Jul 23', 'Aug 23', 'Sept 23', 'Oct 23']

x=sym.symbols('x')
func=0*x
func=sym.Piecewise((func,False),(spline_funcs2[0],x<=X2[1]))

for i in range(1,len(X2)-2):
    func=sym.Piecewise((func,x<=X2[i]),(spline_funcs2[i],x<=X2[i]+1))

func=sym.Piecewise((func,x<=X2[-1]-1),(spline_funcs2[-1],x>X2[-1]-1))

x_new1=np.linspace(12,24,500)
y_new1=[func.subs(x,Y2) for Y2 in x_new1]

plt.plot(x_new1,y_new1)
plt.xticks(X2, x2, rotation=90)
plt.axhline(0, c='black', ls='-') 
plt.scatter(X2, Y2, c='darkgreen') 
plt.title('Cubic Spline Interpolation Post Scandal')
plt.xlabel('Month & Year')
plt.ylabel('Average Views Per Month')
plt.savefig("PostCubicSpline.svg")
plt.show()

#%% linear spline

x,m,b=sym.symbols('x m b')
y=m*x+b
n=len(xs)-1
sx=np.array([])

for i in range(0,n):
    slope=(y-Y[i]).subs(x,xs[i])
    yint=(y-Y[i+1]).subs(x,xs[i+1])
    sol=sym.solve([slope,yint],[m,b])
    sx=np.append(sx,sol[m]*x+sol[b])

print(sx)

xs_new=np.linspace(0,24,500)
ys_new=np.interp(xs_new, xs, Y)
plt.plot(xs_new,ys_new)
plt.xticks(rotation=90)
plt.axhline(0, c='black', ls='-') 
plt.scatter(X, Y, c='darkgreen')
plt.axvline('Oct 22', c='black', ls='--') 
plt.title('Linear Spline Interpolation')
plt.xlabel('Month & Year')
plt.ylabel('Average Views Per Month')
plt.savefig("LinearSpline.svg")
plt.show() 
#%% filled linear spline

x,m,b=sym.symbols('x m b')
y=m*x+b
n=len(xs)-1
sx=np.array([])

for i in range(0,n):
    slope=(y-Y[i]).subs(x,xs[i])
    yint=(y-Y[i+1]).subs(x,xs[i+1])
    sol=sym.solve([slope,yint],[m,b])
    sx=np.append(sx,sol[m]*x+sol[b])

print(sx)

xs_new=np.linspace(0,24,500)
ys_new=np.interp(xs_new, xs, Y)
plt.plot(xs_new,ys_new)
plt.fill_between(xs_new.astype(np.float64),np.asarray(ys_new).astype(np.float64),
                 where=(0<x_new)&(x_new<12),color='forestgreen',alpha=0.2)
plt.fill_between(xs_new.astype(np.float64),np.asarray(ys_new).astype(np.float64),
                 where=(12<x_new)&(x_new<24),color='deeppink',alpha=0.1)
plt.annotate('Area=30,237,110.5',xy=(2,700000),xytext=(2,700000),fontsize=10)
plt.annotate('Area=20,477,398',xy=(15,300000),xytext=(15,300000),fontsize=10)
plt.xticks(rotation=90)
plt.axhline(0, c='black', ls='-') 
plt.scatter(X, Y, c='darkgreen')
plt.axvline('Oct 22', c='black', ls='--') 
plt.title('Linear Spline Interpolation')
plt.xlabel('Month & Year')
plt.ylabel('Average Views Per Month')
plt.savefig("FilledLinearSpline.svg")
plt.show() 

#%% CSR and CTR

def myCSR(f,a,b,n):
    x=sym.symbols('x')
    if n%2!=0:
        print('must be even, try again')
    chop=np.linspace(a,b,n+1)
    h=(b-a)/n
    even=0
    odd=0
    if n%2==0:
        for i in range(2,int(n/2)+1,1):
            even=even+f.subs(x,chop[2*i-2])
        for i in range(1,int(n/2)+1,1):
            odd=odd+f.subs(x,chop[2*i-1])
    csr_sum=(h/3)*(f.subs(x,a)+2*even+4*odd+f.subs(x,b))
    return csr_sum

def myCTR(f,a,b,n):
    x=sym.symbols('x')
    chop=np.linspace(a,b,n+1)
    h=(b-a)/n
    ctr_sum=0
    sol=f.subs(x,a)+f.subs(x,b)
    for i in range(1,n):
        ctr_sum=ctr_sum+f.subs(x,chop[i])
    sol=sol+2*ctr_sum
    sol=(h/2)*sol
    return sol 

#%% split data pre cubic spline CSR

x=sym.symbols('x')
a=11
b=12
print('n Composite Simpsons   Composite Trapezoidal')
for i in range(1,10):
    n=2**i
    CSR=myCSR(spline_funcs1[11], a, b, n)
    print(n,CSR)

#%% split data post cubic spline CSR

x=sym.symbols('x')
a=23
b=24
print('n Composite Simpsons')
for i in range(1,10):
    n=2**i
    CSR=myCSR(spline_funcs2[11], a, b, n)
    print(n,CSR)

#%% cubic spline percent difference

pre_norm=30081881.28/12
post_norm=20179589.12/12
diff=1-(post_norm/pre_norm) 
print(diff)
#%% split data pre linear spline CSR

x=sym.symbols('x')
a=11
b=12
print('n Composite Simpsons')
for i in range(1,2):
    n=2**i
    CSR=myCSR(sx[11], a, b, n)
    print(n,CSR)
    
#%% split data post linear spline CSR

x=sym.symbols('x')
a=23
b=24
print('n Composite Simpsons')
for i in range(1,2):
    n=2**i
    CSR=myCSR(sx[23], a, b, n)
    print(n,CSR)
    
#%% linear spline percent difference

pre_norm=30237110.5/12
post_norm=20477398/12
diff=1-(post_norm/pre_norm) 
print(diff)
#%% linear regression and polynomial regression

xs = np.array([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24]).reshape((-1, 1))
Y = [1764096, 2189474, 3463682, 4465010, 2218769, 2626130, 1505382, 1979853, 
     1636089, 2636005, 2055372, 2439880, 4278833, 2526506, 3611316, 743170, 
     1922541, 1705668, 1326187, 1029268, 1433912, 880781, 1100544, 1191511, 1733155]


from sklearn.linear_model import LinearRegression
lin = LinearRegression()
 
lin.fit(xs, Y)

from sklearn.preprocessing import PolynomialFeatures
 
poly = PolynomialFeatures(degree=24)
X_poly = poly.fit_transform(xs)
 
poly.fit(X_poly, Y)
lin2 = LinearRegression()
lin2.fit(X_poly, Y)

plt.scatter(xs, Y, color='blue')
 
plt.plot(xs, lin.predict(xs), color='red')
plt.title('Linear Regression')
plt.xlabel('Assigned Month')
plt.ylabel('Average Monthly Views')
plt.savefig("LinearRegression.svg")
 
plt.show()

plt.scatter(xs, Y, color='blue') 
 
plt.plot(xs, lin2.predict(poly.fit_transform(xs)),
         color='red')
plt.title('Polynomial Regression')
plt.xlabel('Assigned Month')
plt.ylabel('Average Monthly Views')
plt.savefig("PolynomialRegression.svg")
 
plt.show()

#%% regression split into two parts

X1 = np.array([0,1,2,3,4,5,6,7,8,9,10,11,12]).reshape((-1, 1))
Y1 = [1764096, 2189474, 3463682, 4465010, 2218769, 2626130, 1505382, 1979853, 
     1636089, 2636005, 2055372, 2439880, 4278833]

X2 = np.array([12,13,14,15,16,17,18,19,20,21,22,23,24]).reshape((-1, 1))
Y2 = [4278833, 2526506, 3611316, 743170, 1922541, 1705668, 1326187, 1029268,
      1433912, 880781, 1100544, 1191511, 1733155]

from sklearn.linear_model import LinearRegression
lin = LinearRegression()
 
lin.fit(X1, Y1)

from sklearn.preprocessing import PolynomialFeatures
 
poly = PolynomialFeatures(degree=3)
X_poly = poly.fit_transform(X1)
 
poly.fit(X_poly, Y1)
lin2 = LinearRegression()
lin2.fit(X_poly, Y1)

plt.scatter(X1, Y1, color='blue')
 
plt.plot(X1, lin.predict(X1), color='red')
plt.title('Linear Regression Pre Scandal')
plt.xlabel('Assigned Month')
plt.ylabel('Average Monthly Views')
plt.savefig("PreLinearRegression.svg")
 
plt.show()
print(f"intercept: {lin.intercept_}")
print(f"slope: {lin.coef_}")

plt.scatter(X1, Y1, color='blue')
 
plt.plot(X1, lin2.predict(poly.fit_transform(X1)),
         color='red')
plt.title('Polynomial Regression Pre Scandal')
plt.xlabel('Assigned Month')
plt.ylabel('Average Monthly Views')
plt.savefig("PrePolynomialRegression.svg")
 
plt.show()
print(f"intercept: {lin2.intercept_}")
print(f"slope: {lin2.coef_}")

from sklearn.linear_model import LinearRegression
lin3 = LinearRegression()
 
lin3.fit(X2, Y2)

from sklearn.preprocessing import PolynomialFeatures
 
poly = PolynomialFeatures(degree=2)
X_poly = poly.fit_transform(X2)
 
poly.fit(X_poly, Y2)
lin4 = LinearRegression()
lin4.fit(X_poly, Y2)

plt.scatter(X2, Y2, color='blue')
 
plt.plot(X2, lin3.predict(X2), color='red')
plt.title('Linear Regression Post Scandal')
plt.xlabel('Assigned Month')
plt.ylabel('Average Monthly Views')
plt.savefig("PostLinearRegression.svg")
 
plt.show()
print(f"intercept: {lin3.intercept_}")
print(f"slope: {lin3.coef_}")

plt.scatter(X2, Y2, color='blue')
 
plt.plot(X2, lin4.predict(poly.fit_transform(X2)),
         color='red')
plt.title('Polynomial Regression Post Scandal')
plt.xlabel('Assigned Month')
plt.ylabel('Average Monthly Views')
plt.savefig("PostPolynomialRegression.svg")
 
plt.show()
print(f"intercept: {lin4.intercept_}")
print(f"slope: {lin4.coef_}")
#%% regression split into two parts (filled)

X1 = np.array([0,1,2,3,4,5,6,7,8,9,10,11,12]).reshape((-1, 1))
Y1 = [1764096, 2189474, 3463682, 4465010, 2218769, 2626130, 1505382, 1979853, 
     1636089, 2636005, 2055372, 2439880, 4278833]

X2 = np.array([12,13,14,15,16,17,18,19,20,21,22,23,24]).reshape((-1, 1))
Y2 = [4278833, 2526506, 3611316, 743170, 1922541, 1705668, 1326187, 1029268,
      1433912, 880781, 1100544, 1191511, 1733155]

from sklearn.linear_model import LinearRegression
lin = LinearRegression()
 
lin.fit(X1, Y1)

from sklearn.preprocessing import PolynomialFeatures
 
poly = PolynomialFeatures(degree=3)
X_poly = poly.fit_transform(X1)
 
poly.fit(X_poly, Y1)
lin2 = LinearRegression()
lin2.fit(X_poly, Y1)

plt.scatter(X1, Y1, color='blue')
 
plt.plot(X1, lin.predict(X1), color='red')
plt.title('Linear Regression Pre Scandal')
plt.xlabel('Assigned Month')
plt.ylabel('Average Monthly Views')
plt.fill_between(np.array([0,1,2,3,4,5,6,7,8,9,10,11,12]),lin.predict(X1),
                 0,color='forestgreen',alpha=0.2)
plt.annotate('Area=30,700,223.08',xy=(4,500000),xytext=(4,500000),fontsize=10)
plt.savefig("FilledPreLinearRegression.svg")
 
plt.show()
print(f"intercept: {lin.intercept_}")
print(f"slope: {lin.coef_}")

plt.scatter(X1, Y1, color='blue')
 
plt.plot(X1, lin2.predict(poly.fit_transform(X1)),
         color='red')
plt.title('Polynomial Regression Pre Scandal')
plt.xlabel('Assigned Month')
plt.ylabel('Average Monthly Views')
plt.fill_between(np.array([0,1,2,3,4,5,6,7,8,9,10,11,12]),
                 lin2.predict(poly.fit_transform(X1)),0,color='forestgreen',alpha=0.2)
plt.annotate('Area=30,221,917.49',xy=(4,500000),xytext=(4,500000),fontsize=10)
plt.savefig("FilledPrePolynomialRegression.svg")
 
plt.show()
print(f"intercept: {lin2.intercept_}")
print(f"slope: {lin2.coef_}")

from sklearn.linear_model import LinearRegression
lin3 = LinearRegression()
 
lin3.fit(X2, Y2)

from sklearn.preprocessing import PolynomialFeatures
 
poly = PolynomialFeatures(degree=2)
X_poly = poly.fit_transform(X2)
 
poly.fit(X_poly, Y2)
lin4 = LinearRegression()
lin4.fit(X_poly, Y2)

plt.scatter(X2, Y2, color='blue')
 
plt.plot(X2, lin3.predict(X2), color='red')
plt.title('Linear Regression Post Scandal')
plt.xlabel('Assigned Month')
plt.ylabel('Average Monthly Views')
plt.fill_between(np.array([12,13,14,15,16,17,18,19,20,21,22,23,24]),lin3.predict(X2),
                 0,color='deeppink',alpha=0.1)
plt.annotate('Area=21,676,977.23',xy=(16,500000),xytext=(16,500000),fontsize=10)
plt.savefig("FilledPostLinearRegression.svg")
 
plt.show()
print(f"intercept: {lin3.intercept_}")
print(f"slope: {lin3.coef_}")

plt.scatter(X2, Y2, color='blue')
 
plt.plot(X2, lin4.predict(poly.fit_transform(X2)),
         color='red')
plt.title('Polynomial Regression Post Scandal')
plt.xlabel('Assigned Month')
plt.ylabel('Average Monthly Views')
plt.fill_between(np.array([12,13,14,15,16,17,18,19,20,21,22,23,24]),
                 lin4.predict(poly.fit_transform(X2)),0,color='deeppink',alpha=0.1)
plt.annotate('Area=20,636,650.69',xy=(16,500000),xytext=(16,500000),fontsize=10)
plt.savefig("FilledPostPolynomialRegression.svg")
 
plt.show()
print(f"intercept: {lin4.intercept_}")
print(f"slope: {lin4.coef_}")
#%% split data pre linear regression CSR

x=sym.symbols('x')
f=18728.35164835*x+2445981.813186813
a=0
b=12
print('n Composite Simpsons')
for i in range(1,2):
    n=2**i
    CSR=myCSR(f, a, b, n)
    print(n,CSR)
    
#%% split data post linear regresion CSR

x=sym.symbols('x')
f=-182598.65934066*x+5093190.6373626385
a=12
b=24
print('n Composite Simpsons')
for i in range(1,2):
    n=2**i
    CSR=myCSR(f, a, b, n)
    print(n,CSR)

#%% linear regression percent difference

pre_norm=30700223.08/12
post_norm=21676977.23/12
diff=1-(post_norm/pre_norm) 
print(diff)
#%% split data pre poly regression (cubic) CSR

x=sym.symbols('x')
f=1780431.0082417484+1167936.16978859*x-281160.85389611*x**2+16727.23630536*x**3
a=0
b=12
print('n Composite Simpsons')
for i in range(1,2):
    n=2**i
    CSR=myCSR(f, a, b, n)
    print(n,CSR)
    
#%% split data post poly regresion (quadratic) CSR

x=sym.symbols('x')
f=18530741.746253766-1743088.46553446*x+43346.93906094*x**2
a=12
b=24
print('n Composite Simpsons')
for i in range(1,2):
    n=2**i
    CSR=myCSR(f, a, b, n)
    print(n,CSR)

#%% poly regression percent difference

pre_norm=30221917.49/12
post_norm=20636650.69/12
diff=1-(post_norm/pre_norm) 
print(diff)

#%% logistic model w/ RKMO4 for pre

X1 = np.array([0,1,2,3,4,5,6,7,8,9,10,11,12])
Y1 = [1764096, 2189474, 3463682, 4465010, 2218769, 2626130, 1505382, 1979853, 
     1636089, 2636005, 2055372, 2439880, 4278833]

X2 = np.array([12,13,14,15,16,17,18,19,20,21,22,23,24])
Y2 = [4278833, 2526506, 3611316, 743170, 
     1922541, 1705668, 1326187, 1029268, 1433912, 880781, 1100544, 1191511, 1733155]

# dP/dt=-rP(1-P/K)(1-P/T), P(0)=P_0
# K is carrying capacity (max of data)
# T is threshold population (min of data)
# need to test around different r's

#dP/dt=rp(1-P/K), P(0)=P_0

f1=lambda t,p: -0.15*p*(1-p/4465010)
a=0
b=24
h=0.1
N=int((b-a)/h)
t1=np.linspace(a,b,N+1)
x0=4000000
w1=np.array([x0])

for i in range(0,N):
    k1=h*f1(t1[i],w1[i])
    k2=h*f1(t1[i]+h/2,w1[i]+k1/2)
    k3=h*f1(t1[i]+h/2,w1[i]+k2/2)
    k4=h*f1(t1[i]+h,w1[i]+k3)
    w1=np.append(w1,w1[i]+(1/6)*(k1+2*k2+2*k3+k4))


X1 = np.array([0,1,2,3,4,5,6,7,8,9,10,11,12]).reshape((-1, 1))
Y1 = [1764096, 2189474, 3463682, 4465010, 2218769, 2626130, 1505382, 1979853, 
     1636089, 2636005, 2055372, 2439880, 4278833]

X2 = np.array([13,14,15,16,17,18,19,20,21,22,23,24]).reshape((-1, 1))
Y2 = [2526506, 3611316, 743170, 1922541, 1705668, 1326187, 1029268,
      1433912, 880781, 1100544, 1191511, 1733155]

from sklearn.linear_model import LinearRegression
lin = LinearRegression()
 
lin.fit(X1, Y1)

from sklearn.preprocessing import PolynomialFeatures
 
poly = PolynomialFeatures(degree=3)
X_poly = poly.fit_transform(X1)
 
poly.fit(X_poly, Y1)
lin2 = LinearRegression()
lin2.fit(X_poly, Y1)

lin3 = LinearRegression()
 
lin3.fit(X2, Y2)

from sklearn.preprocessing import PolynomialFeatures
 
poly = PolynomialFeatures(degree=2)
X_poly = poly.fit_transform(X2)
 
poly.fit(X_poly, Y2)
lin4 = LinearRegression()
lin4.fit(X_poly, Y2)

plt.scatter(xs, Y, color='blue')
 
plt.plot(X1, lin.predict(X1), color='red')
#plt.plot(X2, lin3.predict(X2), color='red')
plt.title('Linear Regression vs. Logistic Model')
plt.xlabel('Assigned Month')
plt.ylabel('Average Monthly Views')

 
plt.plot(X2, lin3.predict(X2), color='red')

plt.plot(t1,w1,color='purple')
plt.savefig("RegressionVsLogisticModel.svg")
plt.show()















