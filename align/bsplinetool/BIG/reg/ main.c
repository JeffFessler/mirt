#include		<stdio.h>
#include		<stdlib.h>

#include		"phil.h"
#include		"register.h"
#include		"regFlt3d.h"

int main(void)
{
		struct	rParam		reg;

/* directives.maskCombine
 * During alignment, the masks are geometrically transformed in the same way as the data.
 * The variable 'maskCombine' determines how the mask of the test data and the mask of
 * reference data are to be combined with one another. The registration criterion is
 * computed on the resulting pattern.
 */
		reg.directives.maskCombine = or;
		reg.directives.maskCombine = nor;
		reg.directives.maskCombine = and; /* Strongly recommended choice */
		reg.directives.maskCombine = nand;
		reg.directives.maskCombine = xor;
		reg.directives.maskCombine = nxor;

/* directives.referenceMask
 * The initial mask for the reference data can be specified in three different ways.
 * If 'blank', every voxel is taken into account.
 * If 'provided', the variable 'inMsk2' is a pointer to an arbitrary map of voxels,
 * where the value '0.0' means irrelevant. A relevant voxel is indicated by a non-'0.0'
 * value.
 * If 'computed', the map of relevant voxels is estimated from the data, with the
 * convention that the features of interest are brighter than the background. The
 * algorithm is to perform a very low-pass filtering operation, controlled by the
 * variables 'sx', 'sy' and 'sz'. The result is quantized to two levels (black and white),
 * and the black part is masked out (set to '0.0', irrelevant).
 */
		reg.directives.referenceMask = blank; /* Recommended choice */
		reg.directives.referenceMask = provided;
		reg.directives.referenceMask = computed;

/* directives.testMask
 * The initial mask for the test data is specified in the same way as the mask for the
 * reference data. The same values of the variables 'sx', 'sy' and 'sz' do apply.
 */
		reg.directives.testMask = blank; /* Recommended choice */
		reg.directives.testMask = provided;
		reg.directives.testMask = computed;

/* directives.interpolation
 * The interpolation 'zero' means nearest neighbors. Can be used only when registering
 * the center of gravity of the test and the reference data.
 * 'one' means tri-linear interpolation.
 * 'three' means tri-cubic interpolation.
 */
		reg.directives.interpolation = zero;
		reg.directives.interpolation = one;
		reg.directives.interpolation = three; /* Strongly recommended choice */

/* directives.convergence
 * 'gravity' is a valid option only when neither rotation nor scaling nor skewing nor
 * gray-matching is activated (which leaves only translation). In this mode, the gravity
 * center of the reference is computed (taking the mask into account), and the gravity
 * center of the test is aligned to it. Useful for a gross preliminary alignment when the
 * volumes are offset by more than a dozen voxels. One strategy is to save the
 * displacement resulting from 'gravity' and to use it as first guess for starting the
 * true registration, keeping the test and the reference data unchanged. Another strategy
 * is to substitute the previous test data by the transformed test data resulting from
 * 'gravity', and then to proceed with the true registration. In this case, it is wise to
 * perform 'gravity' with nearest-neighbor interpolation.
 * 'Marquardt' is the true registration algorithm. It is incompatible with nearest-
 * neighbor interpolation.
 */
		reg.directives.convergence = gravity;
		reg.directives.convergence = Marquardt; /* Recommended choice */

/* directives.clipping
 * After the registration has converged, an output volume is created that is the
 * transformed version of the test data supposed to be in alignment with the reference
 * data. Some parts of this output volume might not be relevant, because they might be
 * masked out. If 'clipping' is set to 'TRUE' (or '1'), these irrelevant voxels are
 * rendered with the value given in the variable 'backgrnd'. If 'clipping' is set to
 * 'FALSE' (or '0'), these irrelevant voxels are computed from the data using nearest-
 * neighbor interpolation, irrespective of 'directives.interpolation'.
 */
		reg.directives.clipping = FALSE; /* Recommended choice */
		reg.directives.clipping = TRUE;

/* directives.importFit
 * There are three possibilities for specifying the initial guess of the iterative
 * registration algorithm. If 'importFit' is set to '-1', the file given by the variable
 * 'inFit' is read, and the corresponding transformation is inverted for yielding the
 * initial guess. If 'importFit' is set to '0', the initial guess is the identity
 * transformation, irrespective of the content of the file. If 'importFit' is set to '1',
 * the file yields the initial guess, without inversion.
 */
		reg.directives.importFit = -1;
		reg.directives.importFit = 0; /* Recommended choice */
		reg.directives.importFit = 1;

/* directives.exportFit
 * After the registration has converged, the resulting transformation can be saved in a
 * human-readable file, whose name is given by the variable 'outFit'. If 'exportFit' is
 * set to 'FALSE', nothing happens. If 'exportFit' is set to 'TRUE', the file is
 * (over)written.
 */
		reg.directives.exportFit = FALSE;
		reg.directives.exportFit = TRUE; /* Recommended choice */

/* directives.xTrans, directives.yTrans, directives.zTrans
 * If 'xTrans' is set to 'FALSE', no optimization of the translational parameter along
 * the x-axis is performed, and the value of the initial guess is maintained. If 'xTrans'
 * is set to 'TRUE', optimization takes place. Idem for the y- and z-axis.
 * The recommendation below corresponds to a 3D rigid-body registration task with
 * isotropic voxels.
 */
		reg.directives.xTrans = FALSE;
		reg.directives.xTrans = TRUE; /* Recommended choice */
		reg.directives.yTrans = FALSE;
		reg.directives.yTrans = TRUE; /* Recommended choice */
		reg.directives.zTrans = FALSE;
		reg.directives.zTrans = TRUE; /* Recommended choice */

/* directives.xRot, directives.yRot, directives.zRot
 * If 'xRot' is set to 'FALSE', no optimization of the rotational parameter around
 * the x-axis is performed, and the value of the initial guess is maintained. If 'xRot'
 * is set to 'TRUE', optimization takes place. Idem for the y- and z-axis (an 2D
 * rotation within the x-y plane means a rotation around the z-axis). Note that we
 * consider the voxels to be perfectly cubic. If the reality is different, the result
 * of a pure rotation has no meaning, and it is wiser to consider skewing instead.
 * The recommendation below corresponds to a 3D rigid-body registration task with
 * isotropic voxels.
 */
		reg.directives.xRot = FALSE;
		reg.directives.xRot = TRUE; /* Recommended choice */
		reg.directives.yRot = FALSE;
		reg.directives.yRot = TRUE; /* Recommended choice */
		reg.directives.zRot = FALSE;
		reg.directives.zRot = TRUE; /* Recommended choice */

/* directives.isoScaling
 * If 'isoScaling' is set to 'FALSE', no optimization of the isometric scaling parameter
 * is performed, and the value of the initial guess is maintained. If 'isoScaling' is set
 * to 'TRUE', optimization takes place. Note that we consider the voxels to be perfectly
 * cubic. If the reality is different, the result of an isometric scaling has no meaning,    
 * and it is wiser to consider skewing instead.
 * The recommendation below corresponds to a 3D rigid-body registration task with
 * isotropic voxels.
 */
		reg.directives.isoScaling = FALSE; /* Recommended choice */
		reg.directives.isoScaling = TRUE;

/* directives.xSkew, directives.ySkew, directives.zSkew
 * If 'xSkew' is set to 'FALSE', no optimization of the skewing parameter along the
 * x-axis is performed, and the value of the initial guess is maintained. If 'xSkew'
 * is set to 'TRUE', optimization takes place. Idem for the y- and z-axis. Note that
 * skewing is incompatible with a pure rotation or an isometric scaling.
 * The recommendation below corresponds to a 3D rigid-body registration task with
 * isotropic voxels.
 */
		reg.directives.xSkew = FALSE; /* Recommended choice */
		reg.directives.xSkew = TRUE;
		reg.directives.ySkew = FALSE; /* Recommended choice */
		reg.directives.ySkew = TRUE;
		reg.directives.zSkew = FALSE; /* Recommended choice */
		reg.directives.zSkew = TRUE;

/* directives.matchGrey
 * If 'matchGrey' is set to 'FALSE', no optimization of the gray-level scaling factor
 * is performed, and the value of the initial guess is maintained. If 'matchGrey' is set
 * to 'TRUE', optimization takes place. The optimized parameter is the natural logarithm
 * of the intensity scaling factor; the origin is the value '0.0'; no gray-level offset
 * is optimized.
 */
		reg.directives.matchGrey = FALSE; /* Strongly recommended choice */
		reg.directives.matchGrey = TRUE;

/* directives.greyRendering
 * After the registration has converged, the output volume can be rendered with or without
 * making use of the gray-level scaling parameter, irrespective of its optimization. If
 * 'greyRendering' is set to 'FALSE', no scaling is used (log(1.0) == 0.0), whatever the
 * initial condition was. If 'greyRendering' is set to 'TRUE', the scaling found through
 * optimization (matchGrey == TRUE) or the initial scaling (matchGrey == FALSE) is used for
 * rendering.
 */
		reg.directives.greyRendering = FALSE; /* Recommended choice */
		reg.directives.greyRendering = TRUE;

/* directives.zapMean
 * If 'zapMean' is set to 'FALSE', the input data is left untouched. If zapMean is set
 * to 'TRUE', the test data is modified by removing its average value, and the reference
 * data is also modified by removing its average value prior to optimization.
 */
		reg.directives.zapMean = FALSE; /* Recommended choice */
		reg.directives.zapMean = TRUE;

/* directives.zSqueeze
 * If 'zSqueeze' is set to 'FALSE', the multi-resolution pyramid is computed by performing
 * a size reduction within the x-y plane only, and by leaving the dimension along the
 * z-axis untouched. In this way, volumes of inhomogeneous aspect ratio (in terms of
 * number of pixels per side) can be processed without penalty. If 'zSqueeze' is set to
 * 'TRUE', the size reduction is isotropic. Each additional level brings a reduction by
 * two in each of the three linear dimensions of the volume. If the linear dimension is
 * even, no loss does result; if the linear dimension is odd, the highest coordinate in
 * the volume is first truncated by one before the reduction takes place.
 * The recommendation below corresponds to a 3D rigid-body registration task with
 * isotropic voxels.
 */
		reg.directives.zSqueeze = FALSE;
		reg.directives.zSqueeze = TRUE; /* Recommended choice */

/* nx, ny, nz
 * These variables are used to describe the size of the volumes. All the six volumes
 * (inPtr1, inPtr2, inMsk1, inMsk2, outPtr, mskPtr) must have the same size. In addition,
 * they cannot share memory space. By convention, 2D images require 'nz = 1'.
 * The example below corresponds to a tiny (1 x 1 x 1) volume.
 */
		reg.nx = 1; /* No recommended choice */
		reg.ny = 1; /* No recommended choice */
		reg.nz = 1; /* No recommended choice */

/* inPtr1
 * This pointer should hold the test data. The organization of the data is as follows:
 * contiguous memory, raster order, the fastest running index is x, then y, then z.
 * A typical loop runs thus: {for(;;z++) {for(;;y++) {for(;;x++){} } } }.
 * The test data represents what is aligned. It might be modified by the registration
 * procedure according to the value of the variable 'zapMean'.
 */
		reg.inPtr1 = (float *)malloc(reg.nx * reg.ny * reg.nz * sizeof(float));

/* inPtr2
 * This pointer should hold the reference data: the template toward which the test data
 * is aligned. Might be modified according to 'zapMean'.
 */
		reg.inPtr2 = (float *)malloc(reg.nx * reg.ny * reg.nz * sizeof(float));

/* inMsk1
 * This pointer should hold the test mask if 'testMask == provided', or a mask workspace
 * else. In this latter case, the created mask will be returned.
 */
		reg.inMsk1 = (float *)malloc(reg.nx * reg.ny * reg.nz * sizeof(float));

/* inMsk2
 * This pointer should hold the reference mask if 'referenceMask == provided', or a mask
 * workspace else. In this latter case, the created mask will be returned.
 */
		reg.inMsk2 = (float *)malloc(reg.nx * reg.ny * reg.nz * sizeof(float));

/* outPtr
 * This pointer should provide the space for returning the transformed test data.
 */
		reg.outPtr = (float *)malloc(reg.nx * reg.ny * reg.nz * sizeof(float));

/* mskPtr
 * This pointer should provide the space for returning the transformed test mask. Note
 * that the combination of the reference mask with the transformed mask (as specified by
 * the variable 'directives.maskCombine') is NOT returned.
 */
		reg.mskPtr = (float *)malloc(reg.nx * reg.ny * reg.nz * sizeof(float));

/* inFit
 * Should point to a C-string giving the name of the file from which the initial condition
 * has to be imported, according to the variable 'directives.importFit'.
 */
		reg.inFit = "import.txt";

/* outFit
 * Should point to a C-string giving the name of the file to which the optimal
 * transformation has to be exported, according to the variable 'directives.exportFit'.
 */
		reg.outFit = "export.txt";

/* sx, sy, sz
 * These variables determine the amount of smoothing along the x, y, and z axis,
 * respectively. They give the half-width of a recursive Gaussian smoothing window.
 */
		reg.sx = 4.0; /* No recommended choice */
		reg.sy = 4.0; /* No recommended choice */
		reg.sz = 4.0; /* No recommended choice */

/* firstLambda
 * Marquardt-Levenberg is an adaptative optimizer. The variable 'firstLambda' gives the
 * initial value for the adaptation parameter.
 */
		reg.firstLambda = 1.0; /* Recommended choice */

/* lambdaScale
 * Marquardt-Levenberg is an adaptative optimizer. The variable 'lambdaScale' gives the
 * adaptation step.
 */
		reg.lambdaScale = 4.0; /* Recommended choice */

/* minGain
 * An iterative algorithm needs a convergence criterion. If 'minGain' is set to '0.0',
 * new tries will be performed as long as numerical accuracy permits. If 'minGain'
 * is set between '0.0' and '1.0', the computations will stop earlier, possibly to the
 * price of some loss of accuracy. If 'minGain' is set to '1.0', the algorithm pretends
 * to have reached convergence as early as just after the very first successful attempt.
 */
		reg.minGain = 0.0; /* Strongly recommended choice */

/* epsilon
 * The specification of machine-accuracy is normally machine-dependent. The proposed
 * value has shown good results on a variety of systems; it is the C-constant FLT_EPSILON.
 */
		reg.epsilon = 1.192092896E-07; /* Strongly recommended choice */

/* backgrnd
 * The variable 'backgrnd' holds the arbitrary intensity value to which the masked-out
 * voxels of the transformed test data are rendered when 'directives.clipping' is set to
 * 'TRUE'.
 */
		reg.backgrnd = 0.0; /* No recommended choice */

/* levels
 * This variable specifies how deep the multi-resolution pyramid is. By convention, the
 * finest level is numbered '1', which means that a pyramid of depth '1' is strictly
 * equivalent to no pyramid at all. For best registration results, the rule of thumb is
 * to select a number of levels such that the coarsest representation of the data is a
 * cube between 30 and 60 pixels on each side. To that end, 'directives.zSqueeze' must be
 * judiciously selected.
 */
		reg.levels = 4; /* No recommended choice */

/* lastLevel
 * It is possible to short-cut the optimization before reaching the finest stages, which
 * are the most time-consuming. The variable 'lastLevel' specifies which is the finest
 * level on which optimization is to be performed. If 'lastLevel' is set to the same value
 * as 'levels', the registration will take place on the coarsest stage only. If
 * 'lastLevel' is set to '1', the optimization will take advantage of the whole multi-
 * resolution pyramid.
 */
		reg.lastLevel = 1; /* Strongly recommended choice */

/* Actual call to the procedure */
		regFloat3d(&reg);

		return(0);
}
