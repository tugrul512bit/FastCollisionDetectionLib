#include"FastCollisionDetectionLib.h"
#include"Generator.h"
#include<iostream>

template<typename CoordType>
struct Vector3D
{
	CoordType x,y,z;
	Vector3D<CoordType> crossProduct(Vector3D<CoordType> vec)
	{
		Vector3D<CoordType> res;
		res.x = y*vec.z - z*vec.y;
		res.y = z*vec.x - x*vec.z;
		res.z = x*vec.y - y*vec.x;
		return res;
	}

	Vector3D<CoordType> operator - (Vector3D<CoordType> vec)
	{
		Vector3D<CoordType> result;
		result.x = x-vec.x;
		result.y = y-vec.y;
		result.z = z-vec.z;
		return result;
	}

	Vector3D<CoordType> operator + (Vector3D<CoordType> vec)
	{
		Vector3D<CoordType> result;
		result.x = x+vec.x;
		result.y = y+vec.y;
		result.z = z+vec.z;
		return result;
	}

	Vector3D<CoordType> operator * (CoordType v)
	{
		Vector3D<CoordType> result;
		result.x = x*v;
		result.y = y*v;
		result.z = z*v;
		return result;
	}

	CoordType abs()
	{
		return std::sqrt(x*x+y*y+z*z);
	}

};


template<typename CoordType>
struct PointCloud
{
	CoordType xmin,ymin,zmin;
	CoordType xmax,ymax,zmax;
	Vector3D<CoordType> point[125];
	PointCloud(CoordType x, CoordType y, CoordType z)
	{
		xmin=x-2.5f;
		ymin=y-2.5f;
		zmin=z-2.5f;
		xmax=x-2.5f;
		ymax=y-2.5f;
		zmax=z-2.5f;
		for(int i=0;i<125;i++)
		{
			point[i].x=x+i%5-2.5f;
			point[i].y=y+(i/5)%5-2.5f;
			point[i].z=z+i/25-2.5f;
			if(xmin>point[i].x)
				xmin=point[i].x;
			if(ymin>point[i].y)
				ymin=point[i].y;
			if(zmin>point[i].z)
				zmin=point[i].z;
			if(xmax<point[i].x)
				xmax=point[i].x;
			if(ymax<point[i].y)
				ymax=point[i].y;
			if(zmax<point[i].z)
				zmax=point[i].z;
		}
	}
};

template<typename CoordType>
bool pointCloudIntersection(PointCloud<CoordType>& cl1, PointCloud<CoordType>& cl2)
{
	for(Vector3D<CoordType>& p:cl1.point)
	{
		for(Vector3D<CoordType>& p2:cl2.point)
		{
			if((p-p2).abs()<1.0f)
			{
				return true;
			}
		}
	}
	return false;
}

template<typename CoordType>
bool intersectDim(const CoordType minx, const CoordType maxx, const CoordType minx2, const CoordType maxx2)
{
	return !((maxx < minx2) || (maxx2 < minx));
}

int main()
{

	using cotype = float;
	PointCloud<cotype> ico1(0,0,0);
	// heating the CPU for benchmarking
	for(int i=0;i<10000;i++)
	{
		PointCloud<cotype> ico2(0,0.1f,i*0.1f);
		pointCloudIntersection(ico1,ico2);
	}

	const int N = 10000;
	std::vector<PointCloud<cotype>> objects;
	oofrng::Generator<64> gen;
	for(int i=0;i<N;i++)
	{
		objects.push_back(PointCloud<cotype>(gen.generate1Float()*45,gen.generate1Float()*45,gen.generate1Float()*45));
	}

	// benchmark begin
	size_t nano;
	std::map<int,std::map<int,bool>> collisionMatrix;
	{
		FastColDetLib::Bench bench(&nano);
		for(int i=0;i<N;i++)
			for(int j=0;j<N;j++)
		{
			if(intersectDim(objects[i].xmin,objects[i].xmax,objects[j].xmin,objects[j].xmax))
				if(intersectDim(objects[i].ymin,objects[i].ymax,objects[j].ymin,objects[j].ymax))
					if(intersectDim(objects[i].zmin,objects[i].zmax,objects[j].zmin,objects[j].zmax))
						collisionMatrix[i][j]=pointCloudIntersection(objects[i],objects[j]);
		}
	}
	std::cout<<N*N<<"x collision checks between 2 clouds = "<<nano<<" nanoseconds ("<<(nano/((double)N*N))<<" ns per collision check)"<<std::endl;
	std::cout<<collisionMatrix.size()<<" unique object are in a collision"<<std::endl;
	return 0;

}
