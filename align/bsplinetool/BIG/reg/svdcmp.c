#include		<stddef.h>
#include		<stdlib.h>
#include		<math.h>
#include		"phil.h"

#include		"svdcmp.h"

/************************************************************************/
static	double		pythag			(double				a,
									 double				b);

/************************************************************************/
static	double		pythag			(double				a,
									 double				b) {

		double				absa, absb;

		absa = fabs(a);
		absb = fabs(b);
		if (absa > absb)
		  return(absa * sqrt(1.0 + absb * absb / (absa * absa)));
		else
		  return((absb == 0.0) ? (0.0) : (absb * sqrt(1.0 + absa * absa / (absb * absb))));
} /* End of pythag */

/************************************************************************/
/* FUNCTION: svbksb														*/
/************************************************************************/
extern	int			svbksb			(double				*u,
									 double				w[],
									 double				*v,
									 long				m,
									 long				n,
									 double				b[],
									 double				x[]) {

		double				*tmp;
		double				s;
		long				j, i;

		tmp = (double *)malloc((size_t)n * sizeof(double));
		if (tmp == (double *)NULL) {
		  message("ERROR - Not enough memory for tmp in svbksb");
		  return(ERROR);
		}
		for (j = 0L; (j < n); j++) {
		  s = 0.0;
		  if (w[j] != 0.0) {
			for (i = 0L; (i < m); i++)
			  s += u[i * n + j] * b[i];
			s /= w[j];
		  }
		  tmp[j] = s;
		}
		for (i = 0L; (i < n); i++) {
		  s = 0.0;
		  for (j = 0L; (j < n); j++)
			s += v[i * n + j] * tmp[j];
		  x[i] = s;
		}
		free(tmp);
		return(!ERROR);
} /* End of svbksb */

/************************************************************************/
/* FUNCTION: svdcmp														*/
/************************************************************************/
int					svdcmp			(double				*a,
									 long				m,
									 long				n,
									 double				w[],
									 double				*v) {

		double				*rv1;
		double				anorm, scale;
		double				c, f, g, h, s;
		double				x, y, z;
		long				i, its, j, jj, k, l, nm;
		int					flag;

		rv1 = (double *)malloc((size_t)n * sizeof(double));
		if (rv1 == (double *)NULL) {
		  message("ERROR - Not enough memory for rv1 in svdcmp");
		  return(ERROR);
		}
		g = scale = anorm = 0.0;
		for (i = 0L; (i < n); i++) {
		  l = i + 1L;
		  rv1[i] = scale * g;
		  g = s = scale = 0.0;
		  if (i < m) {
			for (k = i; (k < m); k++)
			  scale += fabs(a[k * n + i]);
			if (scale != 0.0) {
			  for (k = i; (k < m); k++) {
				a[k * n + i] /= scale;
				s += a[k * n + i] * a[k * n + i];
			  }
			  f = a[i * n + i];
			  g = (f >= 0.0) ? (-sqrt(s)) : (sqrt(s));
			  h = f * g - s;
			  a[i * n + i] = f - g;
			  for (j = l; (j < n); j++) {
				for (s = 0.0, k = i; (k < m); k++)
				  s += a[k * n + i] * a[k * n + j];
				f = s / h;
				for (k = i; (k < m); k++)
				  a[k * n + j] += f * a[k * n + i];
			  }
			  for (k = i; (k < m); k++)
				a[k * n + i] *= scale;
			}
		  }
		  w[i] = scale * g;
		  g = s = scale = 0.0;
		  if ((i < m) && (i != (n - 1L))) {
			for (k = l; (k < n); k++)
			  scale += fabs(a[i * n + k]);
			if (scale != 0.0) {
			  for (k = l; (k < n); k++) {
				a[i * n + k] /= scale;
				s += a[i * n + k] * a[i * n + k];
			  }
			  f = a[i * n + l];
			  g = (f >= 0.0) ? (-sqrt(s)) : (sqrt(s));
			  h = f * g - s;
			  a[i * n + l] = f - g;
			  for (k = l; (k < n); k++)
				rv1[k] = a[i * n + k] / h;
			  for (j = l; (j < m); j++) {
				for (s = 0.0, k = l; (k < n); k++)
				  s += a[j * n + k] * a[i * n + k];
				for (k = l; (k < n); k++)
				  a[j * n + k] += s * rv1[k];
			  }
			  for (k = l; (k < n); k++)
				a[i * n + k] *= scale;
			}
		  }
		  anorm = (anorm > (fabs(w[i]) + fabs(rv1[i]))) ? (anorm) : (fabs(w[i]) + fabs(rv1[i]));
		}
		for (i = n - 1L; (i >= 0L); i--) {
		  if (i < (n - 1L)) {
			if (g != 0.0) {
			  for (j = l; (j < n); j++)
				v[j * n + i] = a[i * n + j] / (a[i * n + l] * g);
			  for (j = l; (j < n); j++) {
				for (s = 0.0, k = l; (k < n); k++)
				  s += a[i * n + k] * v[k * n + j];
				for (k = l; (k < n); k++)
				  if (s != 0.0)
					v[k * n + j] += s * v[k * n + i];
			  }
			}
			for (j = l; (j < n); j++)
			  v[i * n + j] = v[j * n + i] = 0.0;
		  }
		  v[i * n + i] = 1.0;
		  g = rv1[i];
		  l = i;
		}
		for (i = (m < n) ? (m - 1L) : (n - 1L); (i >= 0L); i--) {
		  l = i + 1L;
		  g = w[i];
		  for (j = l; (j < n); j++)
			a[i * n + j] = 0.0;
		  if (g != 0.0) {
			g = 1.0 / g;
			for (j = l; (j < n); j++) {
			  for (s = 0.0, k = l; (k < m); k++)
				s += a[k * n + i] * a[k * n + j];
			  f = s * g / a[i * n + i];
			  for (k = i; (k < m); k++)
				if (f != 0.0)
				  a[k * n + j] += f * a[k * n + i];
			}
			for (j = i; (j < m); j++)
			  a[j * n + i] *= g;
		  }
		  else
			for (j = i; (j < m); j++)
			  a[j * n + i] = 0.0;
		  a[i * n + i] += 1.0;
		}
		for (k = n - 1L; (k >= 0L); k--) {
		  for (its = 1L; (its <= (long)SVDiterations); its++) {
			flag = TRUE;
			for ( l = k; (l >= 0L); l--) {
			  nm = l - 1L;
			  if ((fabs(rv1[l]) + anorm) == anorm) {
				flag = FALSE;
				break;
			  }
			  if ((fabs(w[nm]) + anorm) == anorm)
				break;
			}
			if (flag) {
			  c = 0.0;
			  s = 1.0;
			  for (i = l; (i <= k); i++) {
				f = s * rv1[i];
				rv1[i] *= c;
				if ((fabs(f) + anorm) == anorm)
				  break;
				g = w[i];
				h = pythag(f, g);
				w[i] = h;
				h = 1.0 / h;
				c = g * h;
				s = -f * h;
				for (j = 0L; (j < m); j++) {
				  y = a[j * n + nm];
				  z = a[j * n + i];
				  a[j * n + nm] = y * c + z * s;
				  a[j * n + i] = z * c - y * s;
				}
			  }
			}
			z = w[k];
			if (l == k) {
			  if (z < 0.0) {
				w[k] = -z;
				for (j = 0L; (j < n); j++)
				  v[j * n + k] = -v[j * n + k];
			  }
			  break;
			}
			if (its == SVDiterations) {
			  free(rv1);
			  message("ERROR - Singular value decomposition iterations do not converge");
			  return(ERROR);
			}
			x = w[l];
			nm = k - 1L;
			y = w[nm];
			g = rv1[nm];
			h = rv1[k];
			f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
			g = pythag(f, 1.0);
			f = ((x - z) * (x + z) + h * ((y / (f + ((f >= 0.0) ? (fabs(g)) : (-fabs(g)))))
			  - h)) / x;
			c = s = 1.0;
			for (j = l; (j <= nm); j++) {
			  i = j + 1L;
			  g = rv1[i];
			  y = w[i];
			  h = s * g;
			  g = c * g;
			  z = pythag(f, h);
			  rv1[j] = z;
			  c = f / z;
			  s = h / z;
			  f = x * c + g * s;
			  g = g * c - x * s;
			  h = y * s;
			  y *= c;
			  for (jj = 0L; (jj < n); jj++) {
				x = v[jj * n + j];
				z = v[jj * n + i];
				v[jj * n + j] = x * c + z * s;
				v[jj * n + i] = z * c - x * s;
			  }
			  z = pythag(f, h);
			  w[j] = z;
			  if (z != 0.0) {
				z = 1.0 / z;
				c = f * z;
				s = h * z;
			  }
			  f = c * g + s * y;
			  x = c * y - s * g;
			  for (jj = 0L; (jj < m); jj++) {
				y = a[jj * n + j];
				z = a[jj * n + i];
				a[jj * n + j] = y * c + z * s;
				a[jj * n + i] = z * c - y * s;
			  }
			}
			rv1[l] = 0.0;
			rv1[k] = f;
			w[k] = x;
		  }
		}
		free(rv1);
		return(!ERROR);
} /* End of svdcmp */
