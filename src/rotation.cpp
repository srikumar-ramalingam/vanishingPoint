//
// (c) MERL 2012 - 2013
//
/* 
	Created by: Srikumar Ramalingam	
	Rotation estimation from VP
*/

#include "rotation.h"

void rotation::computeRotationFromVP(vp vp_,Eigen::Matrix3d Kmatrix)
{
	
	double vpx_x, vpx_y, vpy_x, vpy_y, vpz_x, vpz_y,norm_1,norm_2,norm_3,lambda_1,lambda_2,lambda_3;
	double rotArray[12][3][3];
	int noRots;
	Eigen::Vector3d vec_1,vec_2,vec_3,tempVec;
	Eigen::Matrix3d tempRot;
	
	vpx_x=vp_.vps[0];
	vpx_y=vp_.vps[1];

	vpy_x=vp_.vps[2];
	vpy_y=vp_.vps[3];

	vpz_x=vp_.vps[4];
	vpz_y=vp_.vps[5];

	// R = inv(K)*[vpx_x vpy_x vpz_x;vpx_y vpy_y vpz_y;1 1 1]*diag(lambda_1,lambda_2,lambda_3]

	// There are four possible solutions and when we fix two vps. 
	// There are 3 possible pairs we can use. 
	// Overall there are 12 possible solutions. 

	// case 1
	noRots=0;
	for(int i=0;i<2;i++)
	for(int j=0;j<2;j++)
	{
		 lambda_1=std::pow(-1.0,i);
		 tempVec(0)=vpx_x;tempVec(1)=vpx_y;tempVec(2)=1;
		 vec_1=Kmatrix.inverse()*tempVec;         
         norm_1=vec_1.norm();

         vec_1=lambda_1*vec_1/norm_1;

         lambda_2=std::pow(-1.0,j);
		 tempVec(0)=vpy_x;tempVec(1)=vpy_y;tempVec(2)=1;
         vec_2=Kmatrix.inverse()*tempVec;
         norm_2=vec_2.norm();
         vec_2=lambda_2*vec_2/norm_2;
		 vec_3=vec_1.cross(vec_2);


		 for(int k=0;k<3;k++)
		 {
			rotArray[noRots][k][0]=vec_1(k);
			rotArray[noRots][k][1]=vec_2(k);
			rotArray[noRots][k][2]=vec_3(k);
		 }
		 noRots++;
	}

	// Case 2 - Rotation computation using VP_Y and VP_Z
	for(int i=0;i<2;i++)
	for(int j=0;j<2;j++)
	{
		 lambda_2=std::pow(-1.0,i);
		 tempVec(0)=vpy_x;tempVec(1)=vpy_y;tempVec(2)=1;
		 vec_2=Kmatrix.inverse()*tempVec;         
         norm_2=vec_2.norm();
         vec_2=lambda_2*vec_2/norm_2;
 
         lambda_3=std::pow(-1.0,j);
		 tempVec(0)=vpz_x;tempVec(1)=vpz_y;tempVec(2)=1;
         vec_3=Kmatrix.inverse()*tempVec;
         norm_3=vec_3.norm();
         vec_3=lambda_3*vec_3/norm_3;
		 vec_1=vec_2.cross(vec_3);

		 for(int k=0;k<3;k++)
		 {
			rotArray[noRots][k][0]=vec_1(k);
			rotArray[noRots][k][1]=vec_2(k);
			rotArray[noRots][k][2]=vec_3(k);
		 }
		 noRots++;
	}
	// Case 2 - Rotation computation using VP_Z and VP_X
	for(int i=0;i<2;i++)
	for(int j=0;j<2;j++)
	{
		 lambda_3=std::pow(-1.0,i);
		 tempVec(0)=vpz_x;tempVec(1)=vpz_y;tempVec(2)=1;
		 vec_3=Kmatrix.inverse()*tempVec;         
         norm_3=vec_3.norm();
         vec_3=lambda_3*vec_3/norm_3;

 
         lambda_1=std::pow(-1.0,j);
		 tempVec(0)=vpx_x;tempVec(1)=vpx_y;tempVec(2)=1;
         vec_1=Kmatrix.inverse()*tempVec;
         norm_1=vec_1.norm();
         vec_1=lambda_1*vec_1/norm_1;
		 vec_2=vec_3.cross(vec_1);

		 for(int k=0;k<3;k++)
		 {
			rotArray[noRots][k][0]=vec_1(k);
			rotArray[noRots][k][1]=vec_2(k);
			rotArray[noRots][k][2]=vec_3(k);

		 }
		 noRots++;
	}
	
	// Find the longest x and y lines. 
	int maxLineIdX,maxLineIdY,lineType;
	double lineLength,maxLengthX,maxLengthY;
	maxLengthX=0;
	maxLengthY=0;
	maxLineIdX=0;
	maxLineIdY=0;
	for(int i=0;i<vp_.noLines;i++)
	{		
		lineType=vp_.getLineType(vp_.p[i*4+0],vp_.p[i*4+1],vp_.p[i*4+2],vp_.p[i*4+3]);

		if ((lineType==0)||(lineType==1))
		{
			lineLength=sqrt((vp_.lines[7*i+0]-vp_.lines[7*i+2])*(vp_.lines[7*i+0]-vp_.lines[7*i+2])+
					 (vp_.lines[7*i+1]-vp_.lines[7*i+3])*(vp_.lines[7*i+1]-vp_.lines[7*i+3]));

			if ((lineType==0)&&(lineLength>maxLengthX))
			{
				maxLengthX=lineLength;
				maxLineIdX=i;
			}

			if ((lineType==1)&&(lineLength>maxLengthY))
			{
				maxLengthY=lineLength;
				maxLineIdY=i;
			}
		}
	}

	// Select one rotation out of all the 12 rotation matrices. 
	Eigen::Vector3d leftEnd,rightEnd,topEnd,bottomEnd;
	if (vp_.lines[7*maxLineIdX]<vp_.lines[7*maxLineIdX+2])
	{
		leftEnd(0)=vp_.lines[7*maxLineIdX];leftEnd(1)=vp_.lines[7*maxLineIdX+1];leftEnd(2)=1;
		rightEnd(0)=vp_.lines[7*maxLineIdX+2];rightEnd(1)=vp_.lines[7*maxLineIdX+3];rightEnd(2)=1;
	}
	else
	{
		leftEnd(0)=vp_.lines[7*maxLineIdX+2];leftEnd(1)=vp_.lines[7*maxLineIdX+3];leftEnd(2)=1;
		rightEnd(0)=vp_.lines[7*maxLineIdX];rightEnd(1)=vp_.lines[7*maxLineIdX+1];rightEnd(2)=1;
	}

	if (vp_.lines[7*maxLineIdX+1]<vp_.lines[7*maxLineIdX+3])
	{
		topEnd(0)=vp_.lines[7*maxLineIdX];topEnd(1)=vp_.lines[7*maxLineIdX+1];topEnd(2)=1;
		bottomEnd(0)=vp_.lines[7*maxLineIdX+2];bottomEnd(1)=vp_.lines[7*maxLineIdX+3];bottomEnd(2)=1;
	}
	else
	{
		topEnd(0)=vp_.lines[7*maxLineIdX+2];topEnd(1)=vp_.lines[7*maxLineIdX+3];topEnd(2)=1;
		bottomEnd(0)=vp_.lines[7*maxLineIdX];bottomEnd(1)=vp_.lines[7*maxLineIdX+1];bottomEnd(2)=1;
	}

	
	int selR=0,check1,check2,check3;
	Eigen::Vector3d leftEndR,rightEndR,topEndR,bottomEndR;

	for(int i=0;i<12;i++)
	{
		for(int j=0;j<3;j++)
			for(int k=0;k<3;k++)
				tempRot(j,k)=rotArray[i][j][k];

		leftEndR=tempRot.transpose()*Kmatrix.inverse()*leftEnd;
		rightEndR=tempRot.transpose()*Kmatrix.inverse()*rightEnd;
		topEndR=tempRot.transpose()*Kmatrix.inverse()*topEnd;
		bottomEndR=tempRot.transpose()*Kmatrix.inverse()*bottomEnd;

		check1=(rightEndR(0)-leftEndR(0))>0;
		check2=(bottomEndR(1)-topEndR(1))>0;    

		check3=(leftEndR(2)>=0)&&(rightEndR(2)>=0)&&(topEndR(2)>=0)&&(bottomEndR(2)>=0);
    
		if (check1&&check2&&check3)
        selR=i;
	}

		printf("Selected rotation:\n");
		for(int j=0;j<3;j++)
		{
			for(int k=0;k<3;k++)
			{
				rot(j,k)=rotArray[selR][j][k];
		//		printf("%lf ",rot(j,k));
			}
		//	printf("\n");
		}

		Eigen::JacobiSVD<Eigen::Matrix3d> svdOfR(rot, Eigen::ComputeFullU | Eigen::ComputeFullV);
		rot = svdOfR.matrixU() * svdOfR.matrixV().transpose();
		for (int j = 0; j<3; j++)
		{
			for (int k = 0; k<3; k++)
			{
				//rot(j, k) = rotArray[selR][j][k];
						printf("%lf ",rot(j,k));
			}
				printf("\n");
		}
}
