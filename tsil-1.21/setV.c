/* Evaluation of V functions. */

#include "internal.h"

/* ******************************************************************* */

TSIL_COMPLEX kUU (TSIL_REAL a,
		  TSIL_REAL b,
		  TSIL_REAL c,
		  TSIL_REAL d,
		  TSIL_COMPLEX s,
		  int vcase)
{
  if (0 == vcase) 
    return (b - a - s)/Delta(s, a, b) + (b - c - d)/Delta(b, c, d) - 1.L/b;

  if ((1 == vcase) || (2 == vcase) || (3 == vcase))
    return 0.L + I*0.L;

  /* DGR commented out in v1.2 */
/*   TSIL_Warn("kUU", "The function V(x,y,z,u) is not implemented in the special case");   */
/*   TSIL_Warn("kUU", "s = (sqrt(x) - sqrt(y)^2) and Delta(y,z,u) = 0.");   */
  return TSIL_Infinity;
}

/* ******************************************************************* */

TSIL_COMPLEX kUT1 (TSIL_REAL a,
		   TSIL_REAL b,
		   TSIL_REAL c,
		   TSIL_REAL d,
		   TSIL_COMPLEX s,
		   int vcase)
{
  if (0 == vcase)
    return 2.L*a*(s - a)/(b*Delta(s, a, b));

  if ((1 == vcase) || (2 == vcase))
    return 4.L*a*(s - a)/(b*Delta(s, a, b));

  if (3 == vcase)
    return -(TSIL_SQRT(a/b) -0.5L)*(
           a*b*b + 2.0L*b*b*(c+d-b) - a*c*c + 2.0L*a*c*d - a*d*d + 
           (b*b - c*c - d*d + 2.0L*c*d)*s)/(s*b*Delta(b,c,d));

  return TSIL_Infinity;
}

/* ******************************************************************* */

TSIL_COMPLEX kUT2 (TSIL_REAL a,
		   TSIL_REAL b,
		   TSIL_REAL c,
		   TSIL_REAL d,
		   TSIL_COMPLEX s,
		   int vcase)
{
  if (0 == vcase)
    return (d/b)*((s - a - b)/Delta(s, a, b) + (b + c - d)/Delta(b, c, d));

  if ((1 == vcase) || (2 == vcase))
    return (-b*b + b*(a + 0.5L*c - 2.5L*d + s) + 2.0L*d*(s - a) + 
            0.5L*(d - c)*(a - s)*(a - s)/b)/(b*Delta(s,a,b));

  if (3 == vcase)
  return d*(b*b + b*(2.5L*s - 0.5L*a - c - d) + 2.0L*s*(c-d)
            + 0.5L*(a-s)*(c-d)*(c-d)/b)/(s*b*Delta(b,c,d));

  return TSIL_Infinity;
}

/* ******************************************************************* */

TSIL_COMPLEX kUS (TSIL_REAL a,
		  TSIL_REAL b,
		  TSIL_REAL c,
		  TSIL_REAL d,
		  TSIL_COMPLEX s,
		  int vcase)
{
  if (0 == vcase) 
    return 2.L*(s - a - b)/(b*Delta(s, a, b));

  if ((1 == vcase) || (2 == vcase))
    return 4.L*(s - a - b)/(b*Delta(s, a, b));

  if (3 == vcase)
    return (2.0L*b*b + b*(s - a - 2.0L*c - 2.0L*d)
            + (a-s)*(c-d)*(c-d)/b)/(s*b*Delta(b,c,d));

  return TSIL_Infinity;
}

/* ******************************************************************* */

TSIL_COMPLEX kUAI (TSIL_REAL a,
		   TSIL_REAL b,
		   TSIL_REAL c,
		   TSIL_REAL d,
		   TSIL_COMPLEX s,
		   int vcase)
{
  if (0 == vcase) 
    return 2.L*(s - a - b)/(b*Delta(s, a, b));

  if ((1 == vcase) || (2 == vcase))
    return 0.0L;

  if (3 == vcase)
    return (b*b + (s-a-b)*(c+d) + (a-s)*(c-d)*(c-d)/b)/(s*b*Delta(b,c,d));

  return TSIL_Infinity;
}

/* ******************************************************************* */

TSIL_COMPLEX kUB (TSIL_REAL a,
		  TSIL_REAL b,
		  TSIL_REAL c,
		  TSIL_REAL d, 
		  TSIL_COMPLEX s,
		  TSIL_REAL QQ,
		  int vcase)
{
  if (0 == vcase)
    return (b*(c -b + d) + (b - c + d)*A(c,QQ) + (b + c - d)*A(d,QQ)
	    )/(b*Delta(b, c, d));

  if ((1 == vcase) || (2 == vcase))
    return -(2.L*a*a - 3.L*a*b + b*b - 4.L*a*s - 3.L*b*s + 2.L*s*s 
           - (A(d,QQ)*((a - b)*(b*b - a*c + a*d) 
           + (b*b + 2.L*a*c + b*c - (2.L*a + b)*d)*s 
           + (d - c)*s*s))/(2.L*b*d) 
           + (A(c,QQ)*((a - b)*(-b*b - a*c + a*d) 
           - (b*b - 2.L*a*c - b*c + (2.L*a + b)*d)*s 
           + (d - c)*s*s))/(2.L*b*c))/(b*Delta(s, a, b));

  if (3 == vcase)
    return 0.0L;

  return TSIL_Infinity;
}

/* ******************************************************************* */

TSIL_COMPLEX kU (TSIL_REAL a,
		 TSIL_REAL b,
		 TSIL_REAL c,
		 TSIL_REAL d, 
		 TSIL_COMPLEX s,
		 TSIL_REAL QQ,
		 int vcase)
{
  TSIL_REAL alphab, alphac, alphad, temp;
  
  if (0 == vcase) 
    return (-s*s/4.L + s*(c + d + 5.L*a/4.L + b/4.L)  
	   - (c + d + a)*(a + b))/(b*Delta(s, a, b))
           + (d + c - b)/Delta(b, c, d);

  if (1 == vcase) {
    if (b > TSIL_TOL) alphab = A(b,QQ)/TSIL_SQRT(b); else alphab = 0.0L;
    if (c > TSIL_TOL) alphac = A(c,QQ)/TSIL_SQRT(c); else alphac = 0.0L;
    if (d > TSIL_TOL) alphad = A(d,QQ)/TSIL_SQRT(d); else alphad = 0.0L;
    return -(4.L*a*a + 3.L*a*b + b*b + 4.L*a*c + 4.L*b*c + 4.L*a*d 
            + 4.L*b*d - 6.5L*a*s - 3.5L*b*s - 4.L*c*s - 4.L*d*s + 2.5L*s*s
            + 2.L*A(a,QQ)*(s - a - 2.L*b + TSIL_SQRT(b)*(alphad + alphac)) 
            + alphab*((s - a - b)*(TSIL_SQRT(b) - alphac - alphad)) 
            + alphac*((b-a-s)*TSIL_SQRT(b) + 4.L*(s-a-b)*TSIL_SQRT(c)) 
            + alphad*((b-a-s)*TSIL_SQRT(b) + 4.L*(s-a-b)*TSIL_SQRT(d)) 
            + 2.0L*alphac*alphad*(s-a-b))/(b*Delta(s,a,b));
  }

  if (2 == vcase) {
    if (c > d) {temp = c; c = d; d = temp;} 
    if (b > TSIL_TOL) alphab = A(b,QQ)/TSIL_SQRT(b); else alphab = 0.0L;
    if (c > TSIL_TOL) alphac = A(c,QQ)/TSIL_SQRT(c); else alphac = 0.0L;
    if (d > TSIL_TOL) alphad = A(d,QQ)/TSIL_SQRT(d); else alphad = 0.0L;
    return -(4.L*a*a + 3.L*a*b + b*b + 4.L*a*c + 4.L*b*c + 4.L*a*d 
            + 4.L*b*d - 6.5L*a*s - 3.5L*b*s - 4.L*c*s - 4.L*d*s + 2.5L*s*s
            + 2.L*A(a,QQ)*(s - a - 2.L*b + TSIL_SQRT(b)*(alphad - alphac)) 
            + alphab*((s - a - b)*(TSIL_SQRT(b) + alphac - alphad)) 
            + alphac*((a-b+s)*TSIL_SQRT(b) + 4.L*(s-a-b)*TSIL_SQRT(c)) 
            - alphad*((a-b+s)*TSIL_SQRT(b) + 4.L*(a+b-s)*TSIL_SQRT(d)) 
            + 2.0L*alphac*alphad*(a+b-s))/(b*Delta(s,a,b));
  }

  if (3 == vcase)
    return -(-0.625L*b*b*b + b*b*(4.L*s - 0.25L*a - c - d)
            + b*(0.375L*a*a + 1.5L*a*c + 1.625L*c*c + 1.5L*a*d 
               + 2.75L*c*d + 1.625L*d*d - 0.375L*a*s - 4.75L*c*s - 4.75L*d*s)
            + (s-a)*(c-d)*(c-d)*(0.375L*a + 0.5L*c + 0.5L*d)/b
            - 0.75L*a*c*c - 0.5L*a*c*d - 0.75L*a*d*d + 0.25L*c*c*s 
            + 1.5L*c*d*s + 0.25L*d*d*s
            + A(a,QQ)*( (c + d - b)*(a*a + (a - 2.0L*b)*(b - s)) 
                        + 2.0L*(b - c + d)*(a - b + s)*A(c,QQ) 
                        + 2.0L*(b + c - d)*(a - b + s)*A(d,QQ))/(2.0L*a)
            + A(b,QQ)*(TSIL_SQRT(a/b) - 1.0L)*(b*(b - c - d) 
                       - (b - c + d)*A(c,QQ) - (b + c - d)*A(d,QQ))
            + (b*b + a*(c + d - b) + (c - 3.L*d)*s - b*(c + d + s))*A(c,QQ) 
            + (b*b + a*(c + d - b) + (d - 3.L*c)*s - b*(c + d + s) + 
              (b + s - a)*A(c,QQ))*A(d,QQ)
            )/(s*b*Delta(b,c,d));

  return TSIL_Infinity;
}

/* ******************************************************************* */

void SetV (TSIL_DATA *foo)
{
  TSIL_REAL    x, y, z, u, v, qq, Ax, Ay, Az, Au, Av;
  TSIL_COMPLEX s, I2zuv, I2xyv, Txuv, Tvxu, Tuxv, Tyzv, Tzyv, Tvyz;
  int vcase;

  x = foo->x;
  y = foo->y;
  z = foo->z;
  u = foo->u;
  v = foo->v;
  s = foo->s;
  qq = foo->qq;


  /* Branch on case for simplicity */
  if (foo->whichFns == STUM) {
    Ax = A(x,qq);
    Ay = A(y,qq);
    Az = A(z,qq);
    Au = A(u,qq);
    Av = A(v,qq);
    I2zuv = I2(z,u,v,qq);
    I2xyv = I2(x,y,v,qq);

    if (x < TSIL_TOL) Txuv = 0.0L + 0.0L*I;
    else Txuv = foo->T[xuv].value;

    if (y < TSIL_TOL) Tyzv = 0.0L + 0.0L*I;
    else Tyzv = foo->T[yzv].value;

    if (z < TSIL_TOL) Tzyv = 0.0L + 0.0L*I;
    else Tzyv = foo->T[zyv].value;

    if (u < TSIL_TOL) Tuxv = 0.0L + 0.0L*I;
    else Tuxv = foo->T[uxv].value;

    if (v < TSIL_TOL) {  
      Tvxu = 0.0L + 0.0L*I;
      Tvyz = 0.0L + 0.0L*I;
    }
    else {
      Tvxu = foo->T[vxu].value;
      Tvyz = foo->T[vyz].value;
    }

    if (z < TSIL_TOL) 
      foo->V[xzuv].value = TSIL_Infinity;
    else if (TSIL_CABS((s - Th2(x,z))/(x*x + z*z)) < TSIL_TOL) 
      foo->V[xzuv].value = TSIL_Infinity;
    else if (0 == Vanalytic(x,z,u,v,s,qq,&(foo->V[xzuv].value)))
      /* This might be an analytic case even if the M integral wasn't! */
      {
	vcase = 0;
	if (TSIL_CABS(1.0L - (u + v + 2.0L*TSIL_SQRT(u*v))/z ) < TSIL_TOL)
	  vcase = 1;
	if (TSIL_CABS(1.0L - (u + v - 2.0L*TSIL_SQRT(u*v))/z ) < TSIL_TOL)
	  vcase = 2;
	if (TSIL_CABS((s-Ps2(x,z))/(x*x + z*z)) < TSIL_TOL)
	  vcase +=3;

	foo->V[xzuv].value = 
	  - kUU(x, z, u, v, s, vcase)*(foo->U[xzuv].value)
	  - kUT1(x, z, u, v, s, vcase)*Txuv
	  - kUT2(x, z, u, v, s, vcase)*Tvxu
	  - kUT2(x, z, v, u, s, vcase)*Tuxv
	  - kUS(x, z, u, v, s, vcase)*(foo->S[uxv].value)
	  + kUAI(x, z, u, v, s, vcase)*(Ax+Au+Av+I2zuv)/2.L
	  - kUB(x, z, u, v, s, qq, vcase)*(foo->B[xz].value)
	  - kU(x, z, u, v, s, qq, vcase);
      }
    
    if (u < TSIL_TOL)
      foo->V[yuzv].value = TSIL_Infinity;
    else if (TSIL_CABS((s - Th2(y,u))/(y*y + u*u)) < TSIL_TOL) 
      foo->V[yuzv].value = TSIL_Infinity;
    else if (0 == Vanalytic(y,u,z,v,s,qq,&(foo->V[yuzv].value)))
      { 
	vcase = 0;
	if (TSIL_CABS(1.0L - (z + v + 2.0L*TSIL_SQRT(z*v))/u ) < TSIL_TOL)
	  vcase = 1;
	if (TSIL_CABS(1.0L - (z + v - 2.0L*TSIL_SQRT(z*v))/u ) < TSIL_TOL)
	  vcase = 2;
	if (TSIL_CABS((s-Ps2(y,u))/(y*y + u*u)) < TSIL_TOL)
	  vcase +=3;
	
	foo->V[yuzv].value =
	  - kUU(y, u, z, v, s, vcase)*(foo->U[yuzv].value)
	  - kUT1(y, u, z, v, s, vcase)*Tyzv
	  - kUT2(y, u, z, v, s, vcase)*Tvyz
	  - kUT2(y, u, v, z, s, vcase)*Tzyv
	  - kUS(y, u, z, v, s, vcase)*(foo->S[vyz].value)
	  + kUAI(y, u, z, v, s, vcase)*(Ay+Az+Av+I2zuv)/2.L
	  - kUB(y, u, z, v, s, qq, vcase)*(foo->B[yu].value)
	  - kU(y, u, z, v, s, qq, vcase);
      }
    
    if (x < TSIL_TOL)
      foo->V[zxyv].value = TSIL_Infinity;
    else if (TSIL_CABS((s - Th2(x,z))/(x*x + z*z)) < TSIL_TOL) 
      foo->V[zxyv].value = TSIL_Infinity;
    else if (0 == Vanalytic(z,x,y,v,s,qq,&(foo->V[zxyv].value)))
      {
	vcase = 0;
	if (TSIL_CABS(1.0L - (y + v + 2.0L*TSIL_SQRT(y*v))/x ) < TSIL_TOL)
	  vcase = 1;
	if (TSIL_CABS(1.0L - (y + v - 2.0L*TSIL_SQRT(y*v))/x ) < TSIL_TOL)
	  vcase = 2;
	if (TSIL_CABS((s-Ps2(x,z))/(x*x + z*z)) < TSIL_TOL)
	  vcase +=3;
	
	foo->V[zxyv].value = 
	  - kUU(z, x, y, v, s, vcase)*(foo->U[zxyv].value)
	  - kUT1(z, x, y, v, s, vcase)*Tzyv
	  - kUT2(z, x, y, v, s, vcase)*Tvyz
	  - kUT2(z, x, v, y, s, vcase)*Tyzv
	  - kUS(z, x, y, v, s, vcase)*(foo->S[vyz].value)
	  + kUAI(z, x, y, v, s, vcase)*(Az+Ay+Av+I2xyv)/2.L
	  - kUB(z, x, y, v, s, qq, vcase)*(foo->B[xz].value)
	  - kU(z, x, y, v, s, qq, vcase);
      }
    
    if (y < TSIL_TOL)
      foo->V[uyxv].value = TSIL_Infinity;
    else if (TSIL_CABS((s - Th2(y,u))/(y*y + u*u)) < TSIL_TOL) 
      foo->V[uyxv].value = TSIL_Infinity;
    else if (0 == Vanalytic(u,y,x,v,s,qq,&(foo->V[uyxv].value)))
      { 
	vcase = 0;
	if (TSIL_CABS(1.0L - (x + v + 2.0L*TSIL_SQRT(x*v))/y ) < TSIL_TOL)
	  vcase = 1;
	if (TSIL_CABS(1.0L - (x + v - 2.0L*TSIL_SQRT(x*v))/y ) < TSIL_TOL)
	  vcase = 2;
	if (TSIL_CABS((s-Ps2(y,u))/(y*y + u*u)) < TSIL_TOL)
	  vcase +=3;
	
	foo->V[uyxv].value = 
	  - kUU(u, y, x, v, s, vcase)*(foo->U[uyxv].value)
	  - kUT1(u, y, x, v, s, vcase)*Tuxv
	  - kUT2(u, y, x, v, s, vcase)*Tvxu
	  - kUT2(u, y, v, x, s, vcase)*Txuv
	  - kUS(u, y, x, v, s, vcase)*(foo->S[uxv].value)
	  + kUAI(u, y, x, v, s, vcase)*(Au+Ax+Av+I2xyv)/2.L
	  - kUB(u, y, x, v, s, qq, vcase)*(foo->B[yu].value)
	  - kU(u, y, x, v, s, qq, vcase);
      }
  }
  else if (foo->whichFns == STU) {

    Ax = A(x,qq);
    Az = A(z,qq);
    Au = A(u,qq);
    Av = A(v,qq);
    I2zuv = I2(z,u,v,qq);

    if (x < TSIL_TOL) Txuv = 0.0L + 0.0L*I;
    else Txuv = foo->T[xuv].value;

    if (u < TSIL_TOL) Tuxv = 0.0L + 0.0L*I;
    else Tuxv = foo->T[uxv].value;

    if (v < TSIL_TOL) {  
      Tvxu = 0.0L + 0.0L*I;
    }
    else {
      Tvxu = foo->T[vxu].value;
    }

    if (z < TSIL_TOL) 
      foo->V[xzuv].value = TSIL_Infinity;
    else if (TSIL_CABS((s - Th2(x,z))/(x*x + z*z)) < TSIL_TOL) 
      foo->V[xzuv].value = TSIL_Infinity;
    else if (0 == Vanalytic(x,z,u,v,s,qq,&(foo->V[xzuv].value)))
      /* This might be an analytic case even if the M integral wasn't! */
      {
	vcase = 0;
	if (TSIL_CABS(1.0L - (u + v + 2.0L*TSIL_SQRT(u*v))/z ) < TSIL_TOL)
	  vcase = 1;
	if (TSIL_CABS(1.0L - (u + v - 2.0L*TSIL_SQRT(u*v))/z ) < TSIL_TOL)
	  vcase = 2;
	if (TSIL_CABS((s-Ps2(x,z))/(x*x + z*z)) < TSIL_TOL)
	  vcase +=3;

	foo->V[xzuv].value = 
	  - kUU(x, z, u, v, s, vcase)*(foo->U[xzuv].value)
	  - kUT1(x, z, u, v, s, vcase)*Txuv
	  - kUT2(x, z, u, v, s, vcase)*Tvxu
	  - kUT2(x, z, v, u, s, vcase)*Tuxv
	  - kUS(x, z, u, v, s, vcase)*(foo->S[uxv].value)
	  + kUAI(x, z, u, v, s, vcase)*(Ax+Au+Av+I2zuv)/2.L
	  - kUB(x, z, u, v, s, qq, vcase)*(foo->B[xz].value)
	  - kU(x, z, u, v, s, qq, vcase);
      }
  }
  else
    ;

  return;
}
