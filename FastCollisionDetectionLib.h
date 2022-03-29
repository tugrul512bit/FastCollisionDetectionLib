/*
 * FastCollisionDetectionLib.h
 *
 *  Created on: Mar 11, 2022
 *      Author: tugrul
 */

#ifndef FASTCOLLISIONDETECTIONLIB_H_
#define FASTCOLLISIONDETECTIONLIB_H_

#include<algorithm>
#include<vector>
#include<map>
#include<unordered_map>
#include<chrono>
#include<memory>
#include<math.h>
#include<queue>
#include<stack>
#include<thread>
#include<mutex>
#include<set>
#include<functional>
#include<condition_variable>
#include<unordered_set>
#include<cmath>
#include<iostream>


namespace FastColDetLib
{

class Bench
{
public:
	Bench(size_t * targetPtr)
	{
		target=targetPtr;
		t1 =  std::chrono::duration_cast< std::chrono::nanoseconds >(std::chrono::high_resolution_clock::now().time_since_epoch());
	}

	~Bench()
	{
		t2 =  std::chrono::duration_cast< std::chrono::nanoseconds >(std::chrono::high_resolution_clock::now().time_since_epoch());
		*target= t2.count() - t1.count();
	}
private:
	size_t * target;
	std::chrono::nanoseconds t1,t2;
};

	template<typename DataType>
	class Memory
	{
	public:
		Memory()
		{
			memory=std::make_shared<std::vector<DataType>>();
			allocPtr=std::make_shared<int>();
			*allocPtr = 0;
			allocPtrPtr=allocPtr.get();
			memory->resize(1024);
			ptr=memory->data();
		}

		inline
		const DataType get(const int index) const noexcept
		{

			return ((DataType* __restrict__ const)ptr)[index];
		}

		inline
		void set(const int index, const DataType data) const noexcept
		{

			((DataType* __restrict__ const)ptr)[index]=data;
		}

		inline
		void readFrom(Memory<DataType>& mem, const int index, const int indexThis, const int n)
		{
			/*
			for(int i=0;i<n;i++)
			{
				ptr[i+indexThis] = mem.ptr[i+index];
			}
			*/

			std::copy(mem.ptr+index,mem.ptr+index+n,ptr+indexThis);
		}

		inline
		const int allocate(const int size)
		{
			const int result = *allocPtrPtr;

			while(size + *allocPtrPtr >= memory->size())
			{
				memory->resize(memory->size()*2);
			}
			*allocPtrPtr += size;
			ptr=memory->data();
			return result;
		}

		inline
		const int capacity()
		{
			return memory->size();
		}

		inline
		void reset()
		{
			*allocPtrPtr = 0;
		}

	private:
		DataType* ptr;
		std::shared_ptr<int> allocPtr;
		int* allocPtrPtr;
		std::shared_ptr<std::vector<DataType>> memory;
	};



	struct MemoryPool
	{
		void clear()
		{
			nodeCollisionMask.reset();
			cellCollisionMask.reset();
			childNodeCount.reset();
			index.reset();
			indexParticle.reset();
			orderParticle.reset();
			minX.reset();
			maxX.reset();
			minY.reset();
			maxY.reset();
			minZ.reset();
			maxZ.reset();
			nodeMinX.reset();
			nodeMinY.reset();
			nodeMinZ.reset();
			nodeInvWidth.reset();
			nodeInvHeight.reset();
			nodeInvDepth.reset();
		}



		// node-particle collision
		Memory<unsigned char> nodeCollisionMask;

		// cell-particle collision
		Memory<unsigned char> cellCollisionMask;
		Memory<char> childNodeCount;
		Memory<int> index;
		Memory<int> indexParticle;
		Memory<int> orderParticle;
		Memory<float> nodeMinX;
		Memory<float> nodeMinY;
		Memory<float> nodeMinZ;
		Memory<float> nodeInvWidth;
		Memory<float> nodeInvHeight;
		Memory<float> nodeInvDepth;
		Memory<float> minX;
		Memory<float> maxX;
		Memory<float> minY;
		Memory<float> maxY;
		Memory<float> minZ;
		Memory<float> maxZ;

		Memory<int> idTmp[8];
		Memory<int> orderTmp[8];


	};

	struct AdaptiveGridV2Fields
	{
		AdaptiveGridV2Fields(MemoryPool mPool, const float minx, const float miny, const float minz,
				const float maxx, const float maxy, const float maxz):mem(mPool),
						minCornerX(minx),minCornerY(miny),minCornerZ(minz),maxCornerX(maxx),maxCornerY(maxy),maxCornerZ(maxz),
						cellWidth    ((maxx-minx)*0.5f),
						cellHeight   ((maxy-miny)*0.5f),
						cellDepth    ((maxz-minz)*0.5f),
						cellWidthInv (1.0f/((maxx-minx)*0.5f)),
						cellHeightInv(1.0f/((maxy-miny)*0.5f)),
						cellDepthInv (1.0f/((maxz-minz)*0.5f))
		{

		}

		MemoryPool mem;
		const float minCornerX;
		const float minCornerY;
		const float minCornerZ;

		const float maxCornerX;
		const float maxCornerY;
		const float maxCornerZ;

		const float cellWidth;
		const float cellHeight;
		const float cellDepth;

		const float cellWidthInv;
		const float cellHeightInv;
		const float cellDepthInv;
	};

	const int testParticleLimit = 32; // maximum particle AABB overlapping allowed on same cell
	class AdaptiveGridV2
	{
	private:


	    // stores a bit in a byte at a position
	    inline void storeBit(unsigned char & data, const unsigned char value, const int pos) noexcept
	    {
	    	data = (value << pos) | (data & ~(1 << pos));
	    }
	public:
		AdaptiveGridV2(MemoryPool mem, const float minx, const float miny, const float minz,
				const float maxx, const float maxy, const float maxz)
		{
			fields = std::make_shared<AdaptiveGridV2Fields>(mem,minx,miny,minz,maxx,maxy,maxz);

		}

		void clear()
		{
			fields->mem.clear();

			// set current (root) node's particle start index to 0
			const int indexParticleStart = fields->mem.index.allocate(1);
			fields->mem.index.set(indexParticleStart,0);

			// set current (root) node's number of particles to 0
			const int indexNumParticles = fields->mem.index.allocate(1);
			fields->mem.index.set(indexNumParticles,0);

			// set current (root) node's child node start
			const int indexChildNodeStart = fields->mem.index.allocate(1);
			fields->mem.index.set(indexChildNodeStart,9);


			// set AABB of current (root) node
			// X
			const int indexBoundMinX = fields->mem.index.allocate(1);
			const int indexBoundMinXFloat = fields->mem.nodeMinX.allocate(1);
			fields->mem.nodeMinX.set(indexBoundMinXFloat,fields->minCornerX);
			fields->mem.index.set(indexBoundMinX,indexBoundMinXFloat);

			// Y
			const int indexBoundMinY = fields->mem.index.allocate(1);
			const int indexBoundMinYFloat = fields->mem.nodeMinY.allocate(1);
			fields->mem.nodeMinY.set(indexBoundMinYFloat,fields->minCornerY);
			fields->mem.index.set(indexBoundMinY,indexBoundMinYFloat);

			// Z
			const int indexBoundMinZ = fields->mem.index.allocate(1);
			const int indexBoundMinZFloat = fields->mem.nodeMinZ.allocate(1);
			fields->mem.nodeMinZ.set(indexBoundMinZFloat,fields->minCornerZ);
			fields->mem.index.set(indexBoundMinZ,indexBoundMinZFloat);

			// cell inverse width
			const int indexWidth = fields->mem.index.allocate(1);
			const int indexWidthFloat = fields->mem.nodeInvWidth.allocate(1);
			fields->mem.nodeInvWidth.set(indexWidthFloat,fields->cellWidthInv);
			fields->mem.index.set(indexWidth,indexWidthFloat);

			// cell inverse height
			const int indexHeight = fields->mem.index.allocate(1);
			const int indexHeightFloat = fields->mem.nodeInvHeight.allocate(1);
			fields->mem.nodeInvHeight.set(indexHeightFloat,fields->cellHeightInv);
			fields->mem.index.set(indexHeight,indexHeightFloat);

			// cell inverse depth
			const int indexDepth = fields->mem.index.allocate(1);
			const int indexDepthFloat = fields->mem.nodeInvDepth.allocate(1);
			fields->mem.nodeInvDepth.set(indexDepthFloat,fields->cellDepthInv);
			fields->mem.index.set(indexDepth,indexDepthFloat);

			fields->mem.childNodeCount.set(fields->mem.childNodeCount.allocate(1),0);
			fields->mem.nodeCollisionMask.set(fields->mem.nodeCollisionMask.allocate(1),0);
		}



		void addParticleDirectAABB(const int id, const float minx, const float miny, const float minz,
				const float maxx, const float maxy, const float maxz,
				int nodeOffset = 0, const bool markStart = false)
		{
			const int particleStart = fields->mem.index.get(nodeOffset);
			const int pId = fields->mem.indexParticle.allocate(1);
			const int oId = fields->mem.orderParticle.allocate(1);
			const int maskId = fields->mem.cellCollisionMask.allocate(1);
			const int maxXId = fields->mem.maxX.allocate(1);
			const int maxYId = fields->mem.maxY.allocate(1);
			const int maxZId = fields->mem.maxZ.allocate(1);
			const int minXId = fields->mem.minX.allocate(1);
			const int minYId = fields->mem.minY.allocate(1);
			const int minZId = fields->mem.minZ.allocate(1);
			fields->mem.indexParticle.set(pId,id);
			fields->mem.orderParticle.set(oId,oId);
			fields->mem.cellCollisionMask.set(maskId,0);
			fields->mem.maxX.set(maxXId,maxx);
			fields->mem.maxY.set(maxYId,maxy);
			fields->mem.maxZ.set(maxZId,maxz);
			fields->mem.minX.set(minXId,minx);
			fields->mem.minY.set(minYId,miny);
			fields->mem.minZ.set(minZId,minz);
			fields->mem.index.set(nodeOffset+1,fields->mem.index.get(nodeOffset+1)+1);
			fields->mem.index.set(nodeOffset+2,0);
			if(pId==0 || markStart)
			{
				fields->mem.index.set(particleStart,pId);

			}
		}




		inline
		const bool intersectDim(const float minx, const float maxx, const float minx2, const float maxx2) const noexcept
		{
			return !((maxx < minx2) || (maxx2 < minx));
		}


		struct NodeTask
		{
			NodeTask(const int n1=0):nodePointer(n1){  }
			const int nodePointer;
		};

		struct LeafTask
		{
			LeafTask(const int n1=0):particlePointer(n1){  }
			int particlePointer;
		};

		// keeps record of unique values inserted
		// works for positive integers (-1 reserved for first comparisons)
		template<typename SignedIntegralType, int n>
		struct FastUnique
		{
		    public:
		    FastUnique()
		    {
		         it=0;
		         for(int i=0;i<n;i++)
		             dict[i]=-1;
		    }

		    inline
		    void insert(const SignedIntegralType val)
		    {
		        const bool result = testImpl(val);
		        dict[it]=(result?val:dict[it]);
		        it+=(result?1:0);
		    }

		    inline
		    const SignedIntegralType get(const int index) const noexcept
		    {
		        return dict[index];
		    }



		    inline
		    const bool test(const SignedIntegralType val) noexcept
		    {
		    	return testImpl(val);
		    }

		    inline
		    const void iterateSet(const SignedIntegralType val) noexcept
		    {
		    	dict[it++]=val;
		    }

		    const int size()
		    {
		        return it;
		    }

		    SignedIntegralType * begin()
		    {
		      return dict;
		    }

		    SignedIntegralType * end()
		    {
		      return dict + it;
		    }

		    private:
		    SignedIntegralType dict[n];
		    SignedIntegralType c[n];
		    int it;

		    inline
		    bool testImpl(const int val) noexcept
		    {
		        for(int i=0;i<n;i++)
		          c[i]=(dict[i]==val);

		        SignedIntegralType s = 0;
		        for(int i=0;i<n;i++)
		          s+=c[i];
		        return s==0;
		    }
		};

		// returns id values of particles
		std::vector<int> findCollisions(const float minx, const float miny, const float minz,
				const float maxx, const float maxy, const float maxz)
		{
			FastUnique<int32_t, testParticleLimit> fastSet;
			std::vector<int> result;
			std::stack<NodeTask> nodesToCompute;
			std::vector<LeafTask> particlesToCompute;


			// push root node to work queue
			nodesToCompute.push(NodeTask(0));



			// traverse all colliding sparse cells
			while(!nodesToCompute.empty() /* stack>=0 */)
			{
				NodeTask task = nodesToCompute.top();
				nodesToCompute.pop();



				const int pointer = fields->mem.index.get(task.nodePointer+2);
				const int numChildNodes = fields->mem.childNodeCount.get(task.nodePointer/9);

				// if this is not a leaf node, traverse all child nodes (they are sparse, so may be less than 8(8bit mask) or 64(64 bit mask))
				if(pointer<0)
				{
					// get current node's information
					const int indexBoundMinX = fields->mem.index.get(task.nodePointer+3);
					const int indexBoundMinY = fields->mem.index.get(task.nodePointer+4);
					const int indexBoundMinZ = fields->mem.index.get(task.nodePointer+5);

					const int indexCellWidthInv = fields->mem.index.get(task.nodePointer+6);
					const int indexCellHeightInv = fields->mem.index.get(task.nodePointer+7);
					const int indexCellDepthInv = fields->mem.index.get(task.nodePointer+8);

					const float minCornerX = fields->mem.nodeMinX.get(indexBoundMinX);
					const float minCornerY = fields->mem.nodeMinY.get(indexBoundMinY);
					const float minCornerZ = fields->mem.nodeMinZ.get(indexBoundMinZ);

					const float cellWidthInv = fields->mem.nodeInvWidth.get(indexCellWidthInv);
					const float cellHeightInv = fields->mem.nodeInvHeight.get(indexCellHeightInv);
					const float cellDepthInv = fields->mem.nodeInvDepth.get(indexCellDepthInv);


					const int indexStartX = std::floor((minx - minCornerX)*cellWidthInv);
					const int indexEndX = std::floor((maxx - minCornerX)*cellWidthInv);

					const int indexStartY = std::floor((miny - minCornerY)*cellHeightInv);
					const int indexEndY = std::floor((maxy - minCornerY)*cellHeightInv);

					const int indexStartZ = std::floor((minz - minCornerZ)*cellDepthInv);
					const int indexEndZ = std::floor((maxz - minCornerZ)*cellDepthInv);


					// prepare cell indicator mask (1 bit = has object, 0 bit = empty))
					unsigned char maskCellsFilled=0;
					for(int k=indexStartZ; k<=indexEndZ; k++)
					{
						if(k<0 || k>=2)
							continue;
						for(int j=indexStartY; j<=indexEndY; j++)
						{
							if(j<0 || j>=2)
								continue;
							for(int i=indexStartX; i<=indexEndX; i++)
							{
								if(i<0 || i>=2)
									continue;

								storeBit(maskCellsFilled,1,i+j*2+k*4);

							}
						}
					}


					const int nodeOffset = -pointer-1;
					for(int i=0;i<numChildNodes;i++)
					{
						// if there is possible collision (accelerated by bit mask for collisions)
						unsigned char cellMask = fields->mem.nodeCollisionMask.get((nodeOffset+i*9)/9);
						if(maskCellsFilled & cellMask)
						{
							nodesToCompute.push(NodeTask(nodeOffset+i*9));
						}
					}
				}
				else
				{
					// this is leaf node

					const int ptr = fields->mem.index.get(task.nodePointer);
					const int n = fields->mem.index.get(task.nodePointer+1);


					const int indexBoundMinX = fields->mem.index.get(task.nodePointer+3);
					const int indexBoundMinY = fields->mem.index.get(task.nodePointer+4);
					const int indexBoundMinZ = fields->mem.index.get(task.nodePointer+5);

					const int indexCellWidthInv = fields->mem.index.get(task.nodePointer+6);
					const int indexCellHeightInv = fields->mem.index.get(task.nodePointer+7);
					const int indexCellDepthInv = fields->mem.index.get(task.nodePointer+8);

					const float minCornerX = fields->mem.nodeMinX.get(indexBoundMinX);
					const float minCornerY = fields->mem.nodeMinY.get(indexBoundMinY);
					const float minCornerZ = fields->mem.nodeMinZ.get(indexBoundMinZ);
					std::cout<<task.nodePointer<<" ? "<<indexBoundMinX<<" ? "<<indexCellWidthInv<<std::endl;
					const float cellWidthInv = fields->mem.nodeInvWidth.get(indexCellWidthInv);
					const float cellHeightInv = fields->mem.nodeInvHeight.get(indexCellHeightInv);
					const float cellDepthInv = fields->mem.nodeInvDepth.get(indexCellDepthInv);

					const int indexStartX = std::floor((minx - minCornerX)*cellWidthInv);
					const int indexEndX = std::floor((maxx - minCornerX)*cellWidthInv);

					const int indexStartY = std::floor((miny - minCornerY)*cellHeightInv);
					const int indexEndY = std::floor((maxy - minCornerY)*cellHeightInv);

					const int indexStartZ = std::floor((minz - minCornerZ)*cellDepthInv);
					const int indexEndZ = std::floor((maxz - minCornerZ)*cellDepthInv);


					// prepare cell indicator mask (1 bit = has object, 0 bit = empty))


					unsigned char maskCellsFilled=0;
					for(int k=indexStartZ; k<=indexEndZ; k++)
					{
						if(k<0 || k>=2)
							continue;
						for(int j=indexStartY; j<=indexEndY; j++)
						{
							if(j<0 || j>=2)
								continue;
							for(int i=indexStartX; i<=indexEndX; i++)
							{
								if(i<0 || i>=2)
									continue;

								storeBit(maskCellsFilled,1,i+j*2+k*4);

							}
						}
					}



					for(int i=0;i<n;i++)
					{
						const int index = ptr+i;
						if(maskCellsFilled & fields->mem.cellCollisionMask.get(index))
						{
							particlesToCompute.push_back(LeafTask(index));
						}

					}



				}
			}

			const int sz = particlesToCompute.size();

			//f(particlesToCompute.data(),sz);
			for(int i=0;i<sz;i++)
			{
				const int index = particlesToCompute[i].particlePointer;

				const int orderId = fields->mem.orderParticle.get(index);
				const int partId = fields->mem.indexParticle.get(orderId);
				if(fastSet.test(partId))
				{



						const float minX = fields->mem.minX.get(orderId);
						const float maxX = fields->mem.maxX.get(orderId);


						if(intersectDim(minx, maxx, minX, maxX))
						{
							const float minY = fields->mem.minY.get(orderId);
							const float maxY = fields->mem.maxY.get(orderId);
							if(intersectDim(miny, maxy, minY, maxY))
							{
								const float minZ = fields->mem.minZ.get(orderId);
								const float maxZ = fields->mem.maxZ.get(orderId);
								if(intersectDim(minz, maxz, minZ, maxZ))
								{
									result.push_back(partId);
									fastSet.iterateSet(partId);
								}
							}
						}

				}
			}

			return result;
		}

		struct PtrN
		{
			int ptr;
			int n;
		};

		// visits all leaf nodes and computes nxn collision pairs
		std::vector<std::pair<int,int>> findCollisionsAll()
		{
			std::unordered_map<int,FastUnique<int32_t, testParticleLimit>> fastSet;
			std::vector<std::pair<int,int>> result;
			std::stack<NodeTask> nodesToCompute;
			std::vector<PtrN> ptrn;


			// push root node to work queue
			nodesToCompute.push(NodeTask(0));



			// traverse all colliding sparse cells
			while(!nodesToCompute.empty() /* stack>=0 */)
			{
				NodeTask task = nodesToCompute.top();
				nodesToCompute.pop();



				const int pointer = fields->mem.index.get(task.nodePointer+2);
				const int numChildNodes = fields->mem.childNodeCount.get(task.nodePointer/9);

				// if this is not a leaf node, traverse all child nodes (they are sparse, so may be less than 8(8bit mask) or 64(64 bit mask))
				if(pointer<0)
				{



					const int nodeOffset = -pointer-1;
					for(int i=0;i<numChildNodes;i++)
					{
						nodesToCompute.push(NodeTask(nodeOffset+i*9));
					}
				}
				else
				{
					// this is leaf node

					const int ptr = fields->mem.index.get(task.nodePointer);
					const int n = fields->mem.index.get(task.nodePointer+1);

					ptrn.push_back({ptr,n});
				}
			}

			// todo: parallelize
			const int sz = ptrn.size();
			for(int k=0;k<sz;k++)
			{
				const int ptr = ptrn[k].ptr;
				const int n = ptrn[k].n;
				for(int i=0;i<n-1;i++)
				{
					const int index = ptr+i;
					const int orderId= fields->mem.orderParticle.get(index);
					const int partId = fields->mem.indexParticle.get(orderId);
					const float minx = fields->mem.minX.get(orderId);
					const float maxx = fields->mem.maxX.get(orderId);
					const float miny = fields->mem.minY.get(orderId);
					const float maxy = fields->mem.maxY.get(orderId);
					const float minz = fields->mem.minZ.get(orderId);
					const float maxz = fields->mem.maxZ.get(orderId);
					auto& map = fastSet[partId];
					for(int j=i;j<n;j++)
					{


						const int index2 = ptr+j;


						const int orderId2 = fields->mem.orderParticle.get(index2);
						const int partId2 = fields->mem.indexParticle.get(orderId2);
						if(partId>=partId2 || !(fields->mem.cellCollisionMask.get(index) & fields->mem.cellCollisionMask.get(index2)))
							continue;
						if(map.test(partId2))
						{

							const float minX = fields->mem.minX.get(orderId2);
							const float maxX = fields->mem.maxX.get(orderId2);

							if(intersectDim(minx, maxx, minX, maxX))
							{

								const float minY = fields->mem.minY.get(orderId2);
								const float maxY = fields->mem.maxY.get(orderId2);
								if(intersectDim(miny, maxy, minY, maxY))
								{

									const float minZ = fields->mem.minZ.get(orderId2);
									const float maxZ = fields->mem.maxZ.get(orderId2);
									if(intersectDim(minz, maxz, minZ, maxZ))
									{
										result.push_back(std::pair<int,int>(partId,partId2));
										map.iterateSet(partId2);
									}
								}
							}
						}
					}
				}
			}
			return result;
		}

		void buildTree()
		{



			int particleStart = fields->mem.index.get(0);
			int numParticle = fields->mem.index.get(1);

			int nodeOffset = 0;



			float minCornerX = fields->mem.nodeMinX.get(fields->mem.index.get(3));
			float minCornerY = fields->mem.nodeMinY.get(fields->mem.index.get(4));
			float minCornerZ = fields->mem.nodeMinZ.get(fields->mem.index.get(5));
			float cellWidthInv =  fields->mem.nodeInvWidth.get(fields->mem.index.get(6));
			float cellHeightInv =  fields->mem.nodeInvHeight.get(fields->mem.index.get(7));
			float cellDepthInv =  fields->mem.nodeInvDepth.get(fields->mem.index.get(8));
			float cellWidth =  1.0f/cellWidthInv;
			float cellHeight =  1.0f/cellHeightInv;
			float cellDepth =  1.0f/cellDepthInv;
			int ctr=0;

			int maxNodeOffset = 9;

			while(nodeOffset <= maxNodeOffset)
			{
				ctr++;

				int ctrTmp[8]={0,0,0,0,0,0,0,0};


				// if child node pointer not set up


				if(fields->mem.index.get(nodeOffset+2)<nodeOffset && fields->mem.index.get(fields->mem.index.get(nodeOffset+2)+2)>=0)
				{
					fields->mem.index.set(fields->mem.index.get(nodeOffset+2)+2,-(nodeOffset+1));
				}


				int childNodeCount = 0;

				if(numParticle > testParticleLimit)
				{


					{

						for(int zz = 0; zz<2; zz++)
							for(int yy = 0; yy<2; yy++)
								for(int xx = 0; xx<2; xx++)
								{
									// allocate node
									const int index0 = xx+yy*2+zz*4;

									fields->mem.orderTmp[index0].reset();
									fields->mem.orderTmp[index0].allocate(numParticle);
								}
					}


					{




						for(int ii=0;ii<numParticle;ii++)
						{

							const int orderParticle = fields->mem.orderParticle.get(particleStart+ii);
							const float minx = fields->mem.minX.get(orderParticle);
							const float miny = fields->mem.minY.get(orderParticle);
							const float minz = fields->mem.minZ.get(orderParticle);
							const float maxx = fields->mem.maxX.get(orderParticle);
							const float maxy = fields->mem.maxY.get(orderParticle);
							const float maxz = fields->mem.maxZ.get(orderParticle);

							const int indexStartX = std::floor((minx - minCornerX)*cellWidthInv);
							const int indexEndX = std::floor((maxx - minCornerX)*cellWidthInv);

							const int indexStartY = std::floor((miny - minCornerY)*cellHeightInv);
							const int indexEndY = std::floor((maxy - minCornerY)*cellHeightInv);

							const int indexStartZ = std::floor((minz -minCornerZ)*cellDepthInv);
							const int indexEndZ = std::floor((maxz - minCornerZ)*cellDepthInv);

							// prepare cell indicator mask (1 bit = has object, 0 bit = empty))

							for(int k=indexStartZ; k<=indexEndZ; k++)
							{
								if(k<0 || k>=2)
									continue;
								for(int j=indexStartY; j<=indexEndY; j++)
								{
									if(j<0 || j>=2)
										continue;
									for(int i=indexStartX; i<=indexEndX; i++)
									{
										if(i<0 || i>=2)
											continue;


										const int index0 = i+j*2+k*4;

										fields->mem.orderTmp[index0].set(ctrTmp[index0],orderParticle);
										ctrTmp[index0]++;
									}
								}
							}
						}


					}


					// add all particles in order (from first child node to last child node)
					childNodeCount=0;
					for(int zz = 0; zz<2; zz++)
						for(int yy = 0; yy<2; yy++)
							for(int xx = 0; xx<2; xx++)
							{
								const int index0 = xx+yy*2+zz*4;
								const int sz = ctrTmp[index0];


								if(sz>0)
								{
									childNodeCount++;

									const int nodeIndexOfs         = fields->mem.index.allocate(9);
									const int particleStartCur     = nodeIndexOfs;
									const int numParticleCur       = nodeIndexOfs+1;
									const int childNodeStartCur    = nodeIndexOfs+2;

									const int nodeBoundMinX        = nodeIndexOfs+3;
									const int nodeBoundMinY        = nodeIndexOfs+4;
									const int nodeBoundMinZ        = nodeIndexOfs+5;

									const int nodeInvWidth         = nodeIndexOfs+6;
									const int nodeInvHeight        = nodeIndexOfs+7;
									const int nodeInvDepth         = nodeIndexOfs+8;
									const int tmpIndex = fields->mem.childNodeCount.allocate(1);
									const int nodeBoundMinXFloat        = fields->mem.nodeMinX.allocate(1);
									const int nodeBoundMinYFloat        = fields->mem.nodeMinY.allocate(1);
									const int nodeBoundMinZFloat        = fields->mem.nodeMinZ.allocate(1);

									const int nodeInvWidthFloat         = fields->mem.nodeInvWidth.allocate(1);
									const int nodeInvHeightFloat        = fields->mem.nodeInvHeight.allocate(1);
									const int nodeInvDepthFloat         = fields->mem.nodeInvDepth.allocate(1);

									fields->mem.nodeMinX.set(nodeBoundMinXFloat,minCornerX+xx*cellWidth);
									fields->mem.nodeMinY.set(nodeBoundMinYFloat,minCornerY+yy*cellHeight);
									fields->mem.nodeMinZ.set(nodeBoundMinZFloat,minCornerZ+zz*cellDepth);

									fields->mem.nodeInvWidth.set(nodeInvWidthFloat,cellWidthInv*2.0f);
									fields->mem.nodeInvHeight.set(nodeInvHeightFloat,cellHeightInv*2.0f);
									fields->mem.nodeInvDepth.set(nodeInvDepthFloat,cellDepthInv*2.0f);


									fields->mem.index.set(nodeBoundMinX,nodeBoundMinXFloat);
									fields->mem.index.set(nodeBoundMinY,nodeBoundMinYFloat);
									fields->mem.index.set(nodeBoundMinZ,nodeBoundMinZFloat);

									fields->mem.index.set(nodeInvWidth,nodeInvWidthFloat);
									fields->mem.index.set(nodeInvHeight,nodeInvHeightFloat);
									fields->mem.index.set(nodeInvDepth,nodeInvDepthFloat);
									const int nodeMaskIndex = fields->mem.nodeCollisionMask.allocate(1);
									unsigned char nodeMask = 0;
									storeBit(nodeMask,1,index0);
									fields->mem.nodeCollisionMask.set(nodeMaskIndex,nodeMask);

									//const int allocOffset = fields->mem.indexParticle.allocate(sz);
									const int allocOffset = fields->mem.orderParticle.allocate(sz);
									fields->mem.cellCollisionMask.allocate(sz);

									//fields->mem.indexParticle.readFrom(fields->mem.idTmp[index0],0,allocOffset,sz);
									fields->mem.orderParticle.readFrom(fields->mem.orderTmp[index0],0,allocOffset,sz);


									fields->mem.index.set(particleStartCur,allocOffset);
									fields->mem.index.set(numParticleCur,sz);
									fields->mem.index.set(childNodeStartCur,nodeOffset);


									maxNodeOffset=particleStartCur;
								}



							}

					fields->mem.childNodeCount.set(nodeOffset/9,childNodeCount);

				}
				else
				{
					fields->mem.childNodeCount.set(nodeOffset/9,0);

					// compute collision mask for each particle

					for(int i=0;i<numParticle;i++)
					{


						const int orderParticle = fields->mem.orderParticle.get(particleStart+i);
						const float minx = fields->mem.minX.get(orderParticle);
						const float miny = fields->mem.minY.get(orderParticle);
						const float minz = fields->mem.minZ.get(orderParticle);
						const float maxx = fields->mem.maxX.get(orderParticle);
						const float maxy = fields->mem.maxY.get(orderParticle);
						const float maxz = fields->mem.maxZ.get(orderParticle);

						const int indexStartX = std::floor((minx - minCornerX)*cellWidthInv);
						const int indexEndX = std::floor((maxx - minCornerX)*cellWidthInv);

						const int indexStartY = std::floor((miny - minCornerY)*cellHeightInv);
						const int indexEndY = std::floor((maxy - minCornerY)*cellHeightInv);

						const int indexStartZ = std::floor((minz -minCornerZ)*cellDepthInv);
						const int indexEndZ = std::floor((maxz - minCornerZ)*cellDepthInv);

						// prepare cell indicator mask (1 bit = has object, 0 bit = empty))
						unsigned char mask = 0;
						for(int kk=indexStartZ; kk<=indexEndZ; kk++)
						{
							if(kk<0 || kk>=2)
								continue;
							for(int jj=indexStartY; jj<=indexEndY; jj++)
							{
								if(jj<0 || jj>=2)
									continue;
								for(int ii=indexStartX; ii<=indexEndX; ii++)
								{
									if(ii<0 || ii>=2)
										continue;

									storeBit(mask,1,ii+jj*2+kk*4);
								}
							}
						}
						fields->mem.cellCollisionMask.set(particleStart+i,mask);
					}
				}

				nodeOffset += 9;
				numParticle=0;
				if(nodeOffset <= maxNodeOffset)
				{
					particleStart = fields->mem.index.get(nodeOffset);
					numParticle = fields->mem.index.get(nodeOffset+1);


					minCornerX = fields->mem.nodeMinX.get(fields->mem.index.get(nodeOffset+3));
					minCornerY = fields->mem.nodeMinY.get(fields->mem.index.get(nodeOffset+4));
					minCornerZ = fields->mem.nodeMinZ.get(fields->mem.index.get(nodeOffset+5));
					cellWidthInv =  fields->mem.nodeInvWidth.get(fields->mem.index.get(nodeOffset+6));
					cellHeightInv =  fields->mem.nodeInvHeight.get(fields->mem.index.get(nodeOffset+7));
					cellDepthInv =  fields->mem.nodeInvDepth.get(fields->mem.index.get(nodeOffset+8));
					cellWidth =  1.0f/cellWidthInv;
					cellHeight =  1.0f/cellHeightInv;
					cellDepth =  1.0f/cellDepthInv;
				}
			}

			std::cout<<"node: "<<ctr<<std::endl;


		}

	private:
		std::shared_ptr<AdaptiveGridV2Fields> fields;
	};





	/*
	 * interface to build various objects that can collide each other
	 *
	 */

	template<typename CoordType>
	class IParticle
	{
	public:
		virtual const CoordType getMaxX()const =0;
		virtual const CoordType getMaxY()const =0;
		virtual const CoordType getMaxZ()const =0;

		virtual const CoordType getMinX()const =0;
		virtual const CoordType getMinY()const =0;
		virtual const CoordType getMinZ()const =0;
		virtual const int getId()const =0;

		const bool intersectX(IParticle<CoordType>* p)
		{
			return !((getMaxX() < p->getMinX()) || (p->getMaxX() < getMinX()));
		}

		const bool intersectY(IParticle<CoordType>* p)
		{
			return !((getMaxY() < p->getMinY()) || (p->getMaxY() < getMinY()));
		}

		const bool intersectZ(IParticle<CoordType>* p)
		{
			return !((getMaxZ() < p->getMinZ()) || (p->getMaxZ() < getMinZ()));
		}

		virtual ~IParticle(){};
	};

	template<typename CoordType>
	class CollisionPair
	{
	public:
		CollisionPair(IParticle<CoordType>* p1Prm=nullptr, IParticle<CoordType>* p2Prm=nullptr)
		{
			p1=p1Prm;
			p2=p2Prm;

		}

		IParticle<CoordType>* getParticle1() const
		{
			return p1;
		}

		IParticle<CoordType>* getParticle2() const
		{
			return p2;
		}
	private:
		IParticle<CoordType> * p1;
		IParticle<CoordType> * p2;
	};


	template<typename CoordType>
	class AdaptiveGrid;


	using GridDataType = char;
	struct MutexWithoutFalseSharing
	{
		std::mutex mut;
		char padding[(64-sizeof(std::mutex))>0?(64-sizeof(std::mutex)):64];
	};

	// Fixed grid of cells (adaptive if a cell overflows)
	template<typename CoordType>
	class FixedGridFields
	{
	public:
		FixedGridFields(const int w, const int h, const int d,	const int s,
				const CoordType minXp, const CoordType minYp, const CoordType minZp,
				const CoordType maxXp, const CoordType maxYp, const CoordType maxZp):
					width(w),height(h),depth(d),widthDiv1(CoordType(1)/w),heightDiv1(CoordType(1)/h),depthDiv1(CoordType(1)/d),storage(s),
					minX(minXp),minY(minYp),minZ(minZp),maxX(maxXp),maxY(maxYp),maxZ(maxZp)
		{


		}

		~FixedGridFields()
		{

		}

		inline
		const int getWidth () const noexcept { return width;};

		inline
		const int getHeight () const noexcept { return height;};

		inline
		const int getDepth () const noexcept { return depth;};

		inline
		const CoordType getWidthDiv1 () const noexcept { return widthDiv1;};

		inline
		const CoordType getHeightDiv1 () const noexcept { return heightDiv1;};

		inline
		const CoordType getDepthDiv1 () const noexcept { return depthDiv1;};

		inline
		const int getStorage () const noexcept { return storage;};


		std::vector<IParticle<CoordType>*> particles;
		std::vector<uint64_t> particlesCollisionMask;

		std::map<IParticle<CoordType>*,std::map<IParticle<CoordType>*,bool>> coll;


		std::map<IParticle<CoordType>*,std::map<IParticle<CoordType>*,bool>> mapping;
		std::vector<CollisionPair<CoordType>> result;
		const int width;
		const int height;
		const int depth;

		const CoordType widthDiv1;
		const CoordType heightDiv1;
		const CoordType depthDiv1;

		const int storage;

		const CoordType minX;
		const CoordType minY;
		const CoordType minZ;
		const CoordType maxX;
		const CoordType maxY;
		const CoordType maxZ;

	};


	template<typename CoordType>
	struct Cmd
	{
		AdaptiveGrid<CoordType>* grid;
		std::mutex* mut;
		std::map<IParticle<CoordType>*,std::map<IParticle<CoordType>*,bool>>* mapping;
		bool* completed;
	};



	template<typename T>
	class SyncQueue
	{
	public:
		SyncQueue(){}
		void push(T t)
		{
			std::unique_lock<std::mutex> lc(m);
			q.push(t);
			c.notify_one();
		}

		void push2(T t)
		{
			std::unique_lock<std::mutex> lc(m);
			q.push(t);
			c.notify_all();
		}

		T pop()
		{
			std::unique_lock<std::mutex> lc(m);
			while(q.empty())
			{
				c.wait(lc);
			}
			T result = q.front();
			q.pop();
			return result;
		}
	private:
		std::queue<T> q;
		std::mutex m;
		std::condition_variable c;
	};


	template<typename CoordType>
	class ThreadPoolFields
	{
	public:
		ThreadPoolFields() {ctr=0; }
		int ctr;
		std::vector<std::thread> worker;
		MutexWithoutFalseSharing mut[7];
		std::vector<int> msg;
		std::vector<std::shared_ptr<SyncQueue<Cmd<CoordType>>>> q;
		~ThreadPoolFields()
		{
			for(unsigned int i=0;i<worker.size();i++)
			{
				std::lock_guard<std::mutex> lg(mut[i].mut);
				msg[i]=0;
				Cmd<CoordType> cmd;
				cmd.grid=nullptr;
				q[i]->push2(cmd);
			}

			for(unsigned int i=0;i<worker.size();i++)
			{

				worker[i].join();
			}
		}
	};

	template<typename CoordType>
	class ThreadPool
	{
	public:
		ThreadPool()
		{
			fields=std::make_shared<ThreadPoolFields<CoordType>>();
			for(int i=0;i<7;i++)
			{
				fields->q.push_back(std::make_shared<SyncQueue<Cmd<CoordType>>>());
				fields->msg.push_back(1);
			}
			auto ptr = fields.get();
			for(int i=0;i<7;i++)
			{

				fields->worker.push_back(std::thread(
						[&,i,ptr]()
						{
							auto fields = ptr;
							bool work = true;
							while(work)
							{
								{
									{
										std::lock_guard<std::mutex> lg(fields->mut[i].mut);
										work=(fields->msg[i]>0);
									}

									Cmd<CoordType> cmd = fields->q[i]->pop();
									if(cmd.grid==nullptr)
										break;
									auto collisions = cmd.grid->getCollisions();

									{
										std::lock_guard<std::mutex> lg(*cmd.mut);
										for(auto& c:collisions)
										{
											(*cmd.mapping)[c.getParticle1()][c.getParticle2()]=true;
										}
										*cmd.completed=true;
									}

								}
							}
						}
				));
			}
		}


		void compute(Cmd<CoordType> cmd)
		{
			fields->q[fields->ctr++%7]->push(cmd);
		}


	private:

		std::shared_ptr<ThreadPoolFields<CoordType>> fields;
	};

	template<typename CoordType>
	class AdaptiveGrid
	{

protected:

		AdaptiveGrid(ThreadPool<CoordType> thr, int depthPrm,
				const CoordType minX, const CoordType minY, const CoordType minZ,
				const CoordType maxX, const CoordType maxY, const CoordType maxZ):thrPool(thr)
		{

			isLeaf = std::make_shared<bool>();
			*isLeaf=true;
			depth = std::make_shared<int>();

			*depth=depthPrm;

			if(*depth<10)
				fields = std::make_shared<FixedGridFields<CoordType>>(4,4,4,300,minX,minY,minZ,maxX,maxY,maxZ);
			else
				fields = std::make_shared<FixedGridFields<CoordType>>(4,4,4,400,minX,minY,minZ,maxX,maxY,maxZ);

			subGrid = std::make_shared<std::vector<AdaptiveGrid<CoordType>>>();

		}



	    // loads a bit from a 8-byte integer at a position
	    inline uint64_t loadBitSizeT(const uint64_t & data, const int pos) noexcept
	    {
	    	return (data>>pos)&1;
	    }

	    // stores a bit in a 8-byte integer at a position
	    inline void storeBitSizeT(uint64_t & data, const uint64_t value, const int pos) noexcept
	    {
	    	data = (value << pos) | (data & ~(((uint64_t)1) << pos));
	    }

public:
		AdaptiveGrid(ThreadPool<CoordType> thr,
				const CoordType minX, const CoordType minY, const CoordType minZ,
				const CoordType maxX, const CoordType maxY, const CoordType maxZ):thrPool(thr)
		{

			isLeaf=std::make_shared<bool>();
			*isLeaf = true;
			depth = std::make_shared<int>();

			*depth=0;

			fields = std::make_shared<FixedGridFields<CoordType>>(4,4,4,100,minX,minY,minZ,maxX,maxY,maxZ);
			subGrid = std::make_shared<std::vector<AdaptiveGrid<CoordType>>>();

		}


		AdaptiveGrid()
		{
			AdaptiveGrid(ThreadPool<CoordType>(),0,0,0,1,1,1);
		}


		void clear()
		{
			*isLeaf=true;
			subGrid->clear();
			fields->particles.clear();
			fields->particlesCollisionMask.clear();
		}


		template<typename Derived>
		void add(Derived * particlesPrm, int n)
		{
			for(int i=0;i<n;i++)
				add(particlesPrm + i);
		}

		// add static particle object pointers to compute all-vs-all comparisons in an optimized way
		// the generated internal data is also used for static vs dynamic collision checking
		template<typename Derived>
		void add(Derived * particlesPrm)
		{

			const int w = fields->getWidth();
			const int h = fields->getHeight();
			const int d = fields->getDepth();
			// grid
			const CoordType xDim = fields->maxX - fields->minX;
			const CoordType yDim = fields->maxY - fields->minY;
			const CoordType zDim = fields->maxZ - fields->minZ;

			// cell
			const CoordType stepX = xDim/w;
			const CoordType stepY = yDim/h;
			const CoordType stepZ = zDim/d;

			const int sto = fields->getStorage();
			const int nPar = fields->particles.size();

			// if current grid leaf is full, convert it to node with 4x4x4 leaves
			if(*isLeaf && (nPar == sto))
			{
				*isLeaf = false;

				// create leaf nodes (4x4x4=64)
				subGrid->reserve(64);
				for(int zz = 0; zz<d; zz++)
					for(int yy = 0; yy<h; yy++)
						for(int xx = 0; xx<w; xx++)
						{

							AdaptiveGrid<CoordType> newGrid(thrPool,*depth+1,fields->minX+stepX*xx,fields->minY+stepY*yy,fields->minZ+stepZ*zz,
									fields->minX+(stepX)*(xx+1),fields->minY+(stepY)*(yy+1),fields->minZ+(stepZ)*(zz+1));



							subGrid->push_back(newGrid);
						}




				for(int ii=0;ii<nPar;ii++)
				{
					// AABB box of particle
					const CoordType minx = fields->particles[ii]->getMinX();
					const CoordType miny = fields->particles[ii]->getMinY();
					const CoordType minz = fields->particles[ii]->getMinZ();

					const CoordType maxx = fields->particles[ii]->getMaxX();
					const CoordType maxy = fields->particles[ii]->getMaxY();
					const CoordType maxz = fields->particles[ii]->getMaxZ();

					const int cellIndexX = std::floor((minx - fields->minX)/stepX);
					const int cellIndexY = std::floor((miny - fields->minY)/stepY);
					const int cellIndexZ = std::floor((minz - fields->minZ)/stepZ);

					const int cellIndexX2 = std::floor((maxx - fields->minX)/stepX);
					const int cellIndexY2 = std::floor((maxy - fields->minY)/stepY);
					const int cellIndexZ2 = std::floor((maxz - fields->minZ)/stepZ);


					for(int zz = cellIndexZ; zz<=cellIndexZ2; zz++)
						for(int yy = cellIndexY; yy<=cellIndexY2; yy++)
							for(int xx = cellIndexX; xx<=cellIndexX2; xx++)
							{
								if(xx<0 || yy<0 || zz<0 || xx>=w || yy>=h || zz>=d)
									continue;

								// overlaps with subgrid, add to it
								(*subGrid)[xx+yy*4+zz*4*4].add(fields->particles[ii]);
							}
				}

				// clear unused particles
				fields->particles.clear();
				fields->particlesCollisionMask.clear();
			}



			{
				// AABB box of particle
				const CoordType minx = (particlesPrm)->getMinX();
				const CoordType miny = (particlesPrm)->getMinY();
				const CoordType minz = (particlesPrm)->getMinZ();

				const CoordType maxx = (particlesPrm)->getMaxX();
				const CoordType maxy = (particlesPrm)->getMaxY();
				const CoordType maxz = (particlesPrm)->getMaxZ();

				const int cellIndexX = std::floor((minx - fields->minX)/stepX);
				const int cellIndexY = std::floor((miny - fields->minY)/stepY);
				const int cellIndexZ = std::floor((minz - fields->minZ)/stepZ);

				const int cellIndexX2 = std::floor((maxx - fields->minX)/stepX);
				const int cellIndexY2 = std::floor((maxy - fields->minY)/stepY);
				const int cellIndexZ2 = std::floor((maxz - fields->minZ)/stepZ);

				uint64_t maskCellsFilled;
				// "gather" operations on neighbor cells should be cache-friendly

				for(int zz = cellIndexZ; zz<=cellIndexZ2; zz++)
					for(int yy = cellIndexY; yy<=cellIndexY2; yy++)
						for(int xx = cellIndexX; xx<=cellIndexX2; xx++)
						{
							if(xx<0 || yy<0 || zz<0 || xx>=w || yy>=h || zz>=d)
								continue;

							storeBitSizeT(maskCellsFilled,1,xx+yy*4+zz*4*4);

							if(!*isLeaf)
							{
								(*subGrid)[xx+yy*4+zz*4*4].add(particlesPrm);
							}
						}

				if(maskCellsFilled)
				{
					if(*isLeaf)
					{
						fields->particlesCollisionMask.push_back(maskCellsFilled);
						fields->particles.push_back(particlesPrm);
					}
				}
			}
		}




		inline
		const bool intersectDim(const CoordType minx, const CoordType maxx, const CoordType minx2, const CoordType maxx2) const noexcept
		{
			return !((maxx < minx2) || (maxx2 < minx));
		}

		// compute collision between given particle and the already-prepared static object grid (after add(..) and getCollisions(..))
		// also returns self-collisions if same particle was added as static particle before (by add(..))
		// thread-safe
		std::vector<IParticle<CoordType>*> getDynamicCollisionListFor(IParticle<CoordType>* particle)
		{
			std::unordered_map<IParticle<CoordType>*,bool> result;
			const int n2 = fields->particles.size();
			result.reserve(n2);

			// AABB box of particle
			const CoordType minx = particle->getMinX();
			const CoordType miny = particle->getMinY();
			const CoordType minz = particle->getMinZ();

			const CoordType maxx = particle->getMaxX();
			const CoordType maxy = particle->getMaxY();
			const CoordType maxz = particle->getMaxZ();

			const CoordType xDim = fields->maxX - fields->minX;
			const CoordType yDim = fields->maxY - fields->minY;
			const CoordType zDim = fields->maxZ - fields->minZ;

			const int w = fields->getWidth();
			const int h = fields->getHeight();
			const int d = fields->getDepth();
			const CoordType stepX = xDim/w;
			const CoordType stepY = yDim/h;
			const CoordType stepZ = zDim/d;


			const int cellIndexX = std::floor((minx - fields->minX)/stepX);
			const int cellIndexY = std::floor((miny - fields->minY)/stepY);
			const int cellIndexZ = std::floor((minz - fields->minZ)/stepZ);

			const int cellIndexX2 = std::floor((maxx - fields->minX)/stepX);
			const int cellIndexY2 = std::floor((maxy - fields->minY)/stepY);
			const int cellIndexZ2 = std::floor((maxz - fields->minZ)/stepZ);
			uint64_t collisionMask=0;
			for(int zz = cellIndexZ; zz<=cellIndexZ2; zz++)
				for(int yy = cellIndexY; yy<=cellIndexY2; yy++)
					for(int xx = cellIndexX; xx<=cellIndexX2; xx++)
					{
						if(xx<0 || yy<0 || zz<0 || xx>=w || yy>=h || zz>=d)
							continue;

						// if selected cell is a cell
						// (if parent is leaf, then it is a cell)
						if(*isLeaf)
						{
							storeBitSizeT(collisionMask,1,xx+yy*4+zz*4*4);
						}
						else // if this is a grid
						{



							auto subResult = (*subGrid)[xx+yy*4+zz*4*4].getDynamicCollisionListFor(particle);

							for(auto& subr:subResult)
							{
								result.emplace(subr,true);
							}
						}
					}

			// if this is a leaf node

			for(int j=0;j<n2;j++)
			{
				if(result.find(fields->particles[j])==result.end())
				if(fields->particlesCollisionMask[j] & collisionMask)
				{

					const CoordType minx2 = fields->particles[j]->getMinX();
					const CoordType maxx2 = fields->particles[j]->getMaxX();
					if(intersectDim(minx, maxx, minx2, maxx2))
					{

						const CoordType miny2 = fields->particles[j]->getMinY();
						const CoordType maxy2 = fields->particles[j]->getMaxY();
						if(intersectDim(miny, maxy, miny2, maxy2))
						{

							const CoordType minz2 = fields->particles[j]->getMinZ();
							const CoordType maxz2 = fields->particles[j]->getMaxZ();
							if(intersectDim(minz, maxz, minz2, maxz2))
							{
								result.emplace(fields->particles[j],true);
							}

						}
					}
				}

			}
			std::vector<IParticle<CoordType>*> resultVec;
			for(auto& res:result)
				resultVec.push_back(res.first);
			return resultVec;
		}

		// returns collision pairs between static objects (and prepares internal data for future dynamic object collision checking), ordered
		std::vector<CollisionPair<CoordType>> getCollisions()
		{
			fields->mapping.clear();
			fields->result.clear();


			const int w = fields->getWidth();
			const int h = fields->getHeight();
			const int d = fields->getDepth();


			// check neighbor cells for a collision of another AABB particle


			std::mutex mut;
			bool completed[64];
			if((!*isLeaf) && *depth>0)
			{
				for(int i=0;i<64;i++)
					completed[i]=false;
			}
			int completedCtr = 0;
			// "gather" operations on neighbor cells should be cache-friendly


			if(!*isLeaf)
			{
				for(int zz = 0; zz<d; zz++)
					for(int yy = 0; yy<h; yy++)
						for(int xx = 0; xx<w; xx++)
				{


						// if at specific layer, enable threads
						if(  (   *(subGrid->data()[xx+yy*4+zz*4*4].isLeaf)   ) && *depth>0)
						{
							Cmd<CoordType> cmd;
							cmd.completed=&completed[completedCtr++];
							cmd.mapping=&fields->mapping;
							cmd.mut=&mut;
							cmd.grid=subGrid->data()+(xx+yy*4+zz*4*4);
							thrPool.compute(cmd);
						}
						else
						{
							auto collisions = (subGrid->data()+(xx+yy*4+zz*4*4))->getCollisions();

							{

								for(auto& c:collisions)
								{
									fields->mapping[c.getParticle1()][c.getParticle2()]=true;
								}

							}
						}

				}
			}


			if(*isLeaf)
			{
				std::map<IParticle<CoordType>*,std::map<IParticle<CoordType>*,bool>> localMap;

				const int nMask = fields->particles.size();
				std::vector<uint64_t> fastTest;

				for(int j=0;j<nMask-1;j++)
				{
					const CoordType minx = fields->particles[j]->getMinX();
					const CoordType maxx = fields->particles[j]->getMaxX();
					const CoordType miny = fields->particles[j]->getMinY();
					const CoordType maxy = fields->particles[j]->getMaxY();
					const CoordType minz = fields->particles[j]->getMinZ();
					const CoordType maxz = fields->particles[j]->getMaxZ();
					for(int i=j+1;i<nMask;i++)
					{
						// if both AABBs in same cell (64bit collision mask = 4x4x4 on/off mapping)
						if(fields->particlesCollisionMask[j] & fields->particlesCollisionMask[i])
						{
							if(fields->particles[j]->getId()<fields->particles[i]->getId())
							{

								const CoordType minx2 = fields->particles[i]->getMinX();
								const CoordType maxx2 = fields->particles[i]->getMaxX();
								if(intersectDim(minx, maxx, minx2, maxx2))
								{

									const CoordType miny2 = fields->particles[i]->getMinY();
									const CoordType maxy2 = fields->particles[i]->getMaxY();
									if(intersectDim(miny, maxy, miny2, maxy2))
									{

										const CoordType minz2 = fields->particles[i]->getMinZ();
										const CoordType maxz2 = fields->particles[i]->getMaxZ();
										if(intersectDim(minz, maxz, minz2, maxz2))
										{
											localMap[fields->particles[j]][fields->particles[i]]=true;
										}
									}
								}
							}
						}
					}
				}


				{
					std::lock_guard<std::mutex> lg(mut);
					for(auto& lm:localMap)
					{
						for(auto& lm2:lm.second)
							fields->mapping[lm.first][lm2.first]=true;
					}
				}
			}

			// if at specific layer, wait for threads
			if((!*isLeaf) && *depth>0)
			{
				bool comp = false;
				while(!comp)
				{
					comp=true;
					{
						std::lock_guard<std::mutex> lg(mut);
						for(int cmd=0;cmd<completedCtr;cmd++)
						{
							comp&=completed[cmd];
						}
					}
					std::this_thread::yield();
				}
			}

			for(auto& m:fields->mapping)
			{
				for(auto& m2:m.second)
				{

					fields->result.push_back(CollisionPair<CoordType>(m.first,m2.first));
				}
			}




			return fields->result;
		}
private:
		std::shared_ptr<FixedGridFields<CoordType>> fields;
		std::shared_ptr<std::vector<AdaptiveGrid<CoordType>>> subGrid;
		ThreadPool<CoordType> thrPool;
		std::shared_ptr<int> depth;
		std::shared_ptr<bool> isLeaf;

	};


	// axis-aligned bounding-box collision detection
	template<typename CoordType>
	class BruteForce
	{
public:
		BruteForce()
		{

		}

		template<typename Derived>
		void add(Derived * particlesPrm, int numParticlesToAdd)
		{
			for(int i=0;i<numParticlesToAdd;i++)
				particles.push_back(static_cast<IParticle<CoordType>*>(particlesPrm+i));
		}

		std::vector<CollisionPair<CoordType>> getCollisions()
		{
			std::vector<CollisionPair<CoordType>> result;
			idMap.clear();
			collisionPairs.clear();
			const int sz = particles.size();
			for(int i=0;i<sz-1;i++)
			{
				idMap[particles[i]->getId()]=particles[i];

				for(int j=i+1;j<sz;j++)
				{

					if( particles[i]->intersectX(particles[j]) && particles[i]->intersectY(particles[j]) && particles[i]->intersectZ(particles[j]))
					{

						collisionPairs.push_back(CollisionPair<CoordType>(particles[i],particles[j]));
					}
				}
			}
			std::sort(collisionPairs.begin(),collisionPairs.end(),[](CollisionPair<CoordType>& c1, CollisionPair<CoordType>& c2){
				return c1.getParticle1()->getId()<c2.getParticle1()->getId();
			});
			result=collisionPairs;
			return result;
		}
private:
		std::vector<IParticle<CoordType>*> particles;
		std::vector<CollisionPair<CoordType>> collisionPairs;
		std::map<int,IParticle<CoordType>*> idMap;
	};
}



#endif /* FASTCOLLISIONDETECTIONLIB_H_ */
