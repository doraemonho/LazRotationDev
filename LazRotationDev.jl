module LazRotationDev
using PyCall, Images, HDF5, Statistics, StatsBase;
using LazCore, LazType

##############################################################################
#
# Copyright (c) 2020
# Ka Wai Ho, Ka Ho Yuen, Yue Hu, Junda Chen and Alex Lazarian
# All Rights Reserved.
#
# ​This program is free software: you can redistribute it and/or modify
# ​it under the terms of the GNU General Public License as published by
# ​the Free Software Foundation, either version 3 of the License, or
# ​(at your option) any later version.
# ​
# This program is distributed in the hope that it will be useful,
# ​but WITHOUT ANY WARRANTY; without even the implied warranty of
# ​MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# ​GNU General Public License for more details.
# ​You should have received a copy of the GNU General Public License
# ​along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
##############################################################################

###############################################################################
#                                                                             #
#                                    Code for                                 # 
#                                                                             #
#                             Peroidic Cube Rotation                          #
#                                                                             #
###############################################################################
#   
#
#       Version     : v1.0.1(10/11/2020)
#       Author      : KW HO @ Lazarian Technology
#       Description : Basic Vector and Scalar Rotation with linear interpolation 
#                     in MHD simulation. 
#       Capability  : Work for julia version 1.2.0
#       Cuation     : For periodic Cubes only!
#
###############################################################################


export VectorRotation
export ScalarRotation


function RotationMatrixGenerator(θx::Number,θy::Number,θz::Number)
	Tx = [      1        0       0
		        0  cos(θx) -sin(θx)
		        0  sin(θx) +cos(θx)];
	Ty = [cos(θy)        0 +sin(θy)
	     		0        1       0	  
		  -sin(θy)       0 +cos(θy)];
	Tz = [cos(θz) -sin(θz)       0
		  sin(θz) +cos(θz)       0
		       0         0       1 ];

	return inv(Tx*Ty*Tz)
end

function VolumeScalerIntp(x_fl::Number,y_fl::Number,z_fl::Number,
						  x_cl::Number,y_cl::Number,z_cl::Number,
						  Δx::Number,Δy::Number,Δz::Number,d::Cube)

	c00 = (1-Δx)*d[x_fl,y_fl,z_fl] + Δx*d[x_cl,y_fl,z_fl];
	c01 = (1-Δx)*d[x_fl,y_fl,z_cl] + Δx*d[x_cl,y_fl,z_cl];
	c10 = (1-Δx)*d[x_fl,y_cl,z_fl] + Δx*d[x_cl,y_cl,z_fl];
	c11 = (1-Δx)*d[x_cl,y_fl,z_fl] + Δx*d[x_cl,y_cl,z_cl];

	c0  = c00*(1-Δy)+c10*Δy;
	c1  = c10*(1-Δy)+c11*Δy;

	c   = c0*(1-Δz)+c1*Δz;

	return c

end

function PeroidericChecking(X::Vec,N::Vec)
	# Peroideric Checking
	for k = 1:3
		if X[k] > N[k]
			X[k] = X[k] - N[k] + 1;
		elseif X[k] < 1
			X[k] = X[k] + N[k] - 1;
		end
	end
end


function ScalarRotation(d::Cube,θx::Number,θy::Number,θz::Number;backward_compatible=false)
    if backward_compatible == true
        θx,θy,θz = -θx,-θy,-θz;
    end
    
    #Nx=Ny=Nz = even number version
	nx,ny,nz = size(d);
	R = zeros(size(d));

	T = RotationMatrixGenerator(θx,θy,θz);
	Y = [0.0, 0.0, 0.0];
	

	N      = [nx,ny,nz];
    NHalf  = [nx/2,ny/2,nz/2].+0.5;
    
	@inbounds for k = 1:nz
		Y[3] = k-NHalf[3];
		for j = 1:ny
			Y[2] = j-NHalf[2];
			for i = 1:nx
				Y[1] = i-NHalf[1];
				X    = T*Y;
				X[1],X[2],X[3] = X[1]+NHalf[1],X[2]+NHalf[2],X[3]+NHalf[3];\
                PeroidericChecking(X,N);
                x_fl,y_fl,z_fl = Int(floor(X[1])),Int(floor(X[2])),Int(floor(X[3]));
                x_cl,y_cl,z_cl =  Int(ceil(X[1])), Int(ceil(X[2])), Int(ceil(X[3]));

                Δx = X[1]-x_fl;
                Δy = X[2]-y_fl;
                Δz = X[3]-z_fl;

                #linear intp;
                R[i,j,k] = VolumeScalerIntp(x_fl,y_fl,z_fl,
                                            x_cl,y_cl,z_cl,
                                    	    Δx,Δy,Δz,d);
			end
		end
	end
    return R
end

function VectorRotation(iv::Cube,jv::Cube,kv::Cube,θx::Number,θy::Number,θz::Number;backward_compatible=false)
    if backward_compatible == true
        θx,θy,θz = -θx,-θy,-θz;
    end
    
    #Nx=Ny=Nz = even number version
	nx,ny,nz = size(iv);
	Riv = zeros(size(iv));
	Rjv = zeros(size(iv));
	Rkv = zeros(size(iv));

	T = RotationMatrixGenerator(θx,θy,θz);
	Y = [0.0, 0.0, 0.0];
	

	N      = [nx,ny,nz];
    NHalf  = [nx/2,ny/2,nz/2].+0.5;
    VI     = zeros(size(N));
    RVI    = zeros(size(N));
	@inbounds for k = 1:nz
		Y[3] = k-NHalf[3];
		for j = 1:ny
			Y[2] = j-NHalf[2];
			for i = 1:nx
				Y[1] = i-NHalf[1];
				X    = T*Y;
				X[1],X[2],X[3] = X[1]+NHalf[1],X[2]+NHalf[2],X[3]+NHalf[3];
                    PeroidericChecking(X,N);
                    x_fl,y_fl,z_fl = Int(floor(X[1])),Int(floor(X[2])),Int(floor(X[3]));
                    x_cl,y_cl,z_cl =  Int(ceil(X[1])), Int(ceil(X[2])), Int(ceil(X[3]));

                    Δx = X[1]-x_fl;
                    Δy = X[2]-y_fl;
                    Δz = X[3]-z_fl;

                    #linear intp;
                    VI[1] = VolumeScalerIntp(x_fl,y_fl,z_fl,
                                                x_cl,y_cl,z_cl,
                                                Δx,Δy,Δz,iv);
                    VI[2] = VolumeScalerIntp(x_fl,y_fl,z_fl,
					                            x_cl,y_cl,z_cl,
					                            Δx,Δy,Δz,jv);
                    VI[3] = VolumeScalerIntp(x_fl,y_fl,z_fl,
					                            x_cl,y_cl,z_cl,
					                            Δx,Δy,Δz,kv);
                    RVI = T*VI;

                    Riv[i,j,k] = RVI[1];
                    Rjv[i,j,k] = RVI[2];
                    Rkv[i,j,k] = RVI[3];

			end
		end
	end

	return Riv,Rjv,Rkv

end



end
