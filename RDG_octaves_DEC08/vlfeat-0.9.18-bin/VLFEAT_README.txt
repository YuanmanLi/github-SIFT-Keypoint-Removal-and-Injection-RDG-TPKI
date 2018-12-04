%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
The corresponding version of Version is vlfeat-0.9.18. 
More information please refer to http://www.vlfeat.org/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This is recompiled using VS2010, other version may not be supported
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Here only the sift component is used. 
In order to get some intermediate results of vl_sift, which are not provided by the offical version,
we have made the following changes. However, our changes do not affect anything of the genuine one. 
One can direcly get the same frames and descriptors, and also some other additional information.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Input paramters:
'outputOctave':  To get all the keypoints in a given octave. Note that those keypoints are without interpolation. (optional)

'ThreEqual':  the threshold to justify whether two DoG vaules are equal. They are equal only when the differece is smaller than ThreEqual (default: 0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Output parameters: [frame, descriptor, keypoints, keyrep]
keypoints: 3-dimensional vector, which is useful only when outputOctave parameter is given. (dimension 3)
keyrep:    1-dimensional vector, the map of the duplicated keypoints. if one keypoint is not an duplicated one, then assign value 1.    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Interested readers can send us emails for the code.   (yuanmanx.li@gmail.com)
