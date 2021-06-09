Hi Sylvia,

Thanks for your prompt reply and useful comments.
We added the citations you suggested: a sentence "Other studies have been done in this direction, mostly using a parametric model for the PSD" to acknowledge other works.

Regarding the likelihood term in equation 26:
	-I don't see a problem in the lhs being TD and the rhs being FD: the likelihood in TD can depend on its FD version
	- We fixed the square: instead of $(\cdot)^2$, we wrote $|\cdot|^2$
	- The integral over the PSD is the normalization of the gaussian. In a finite dimensional case, it corresponds to the determinant of the covariance matrix A, appearing as $1/\sqrt{det(A)}$. This also explains the factor 0.5. This term does not matter if the PSD does not depend on $\theta$ (it cancels out in the normalization), but it does get important if it depends on $\theta$.

We uploaded the changes in the v3 version of the file on DCC

Best,

Stefano
