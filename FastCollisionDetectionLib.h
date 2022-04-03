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



	/*
	 * interface to build various objects that can collide each other
	 *
	 */


	inline
	const int intersectDim(const float minx, const float maxx, const float minx2, const float maxx2) noexcept
	{
		return !((maxx < minx2) || (maxx2 < minx));
	}

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
	    void reset()
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
		DataType * getPtr(const int index) const noexcept
		{

			return ptr+index;
		}

		inline
		DataType& getRef(const int index) const noexcept
		{

			return ((DataType* __restrict__ const)ptr)[index];
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
			std::copy(mem.ptr+index,mem.ptr+index+n,ptr+indexThis);
		}

		inline
		void writeTo(std::vector<DataType>& vec)
		{
			std::copy(ptr,ptr+*allocPtrPtr,vec.data());
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
		const int size()
		{
			return *allocPtrPtr;
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


	constexpr int testParticleLimit = 32; // maximum particle AABB overlapping allowed on same cell
	struct MemoryPool
	{
		void clear()
		{
			nodeCollisionMask.reset();

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
			leafOffset.reset();
		}



		// node-particle collision
		Memory<uint64_t> nodeCollisionMask;

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

		Memory<int> idTmp[64];
		Memory<int> orderTmp[64];

		Memory<std::pair<int,int>> allPairsColl;

		Memory<FastUnique<int32_t, testParticleLimit>> allPairsCollmapping;
		Memory<int> leafOffset;

	};

	struct AdaptiveGridV2Fields
	{
		AdaptiveGridV2Fields(MemoryPool mPool, const float minx, const float miny, const float minz,
				const float maxx, const float maxy, const float maxz):mem(mPool),
						minCornerX(minx),minCornerY(miny),minCornerZ(minz),maxCornerX(maxx),maxCornerY(maxy),maxCornerZ(maxz),
						cellWidth    ((maxx-minx)*0.25f),
						cellHeight   ((maxy-miny)*0.25f),
						cellDepth    ((maxz-minz)*0.25f),
						cellWidthInv (1.0f/((maxx-minx)*0.25f)),
						cellHeightInv(1.0f/((maxy-miny)*0.25f)),
						cellDepthInv (1.0f/((maxz-minz)*0.25f))

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


	class AdaptiveGridV2
	{
	private:


	    // stores a bit in a byte at a position
	    inline void storeBit(uint64_t & data, const uint64_t value, const int pos) noexcept
	    {
	    	data = (value << pos) | (data & ~(((uint64_t)1) << pos));
	    }



		void comp4vs4(	const int * const __restrict__ partId1, const int * const __restrict__ partId2,
							const float * const __restrict__ minx1, const float * const __restrict__ minx2,
							const float * const __restrict__ miny1, const float * const __restrict__ miny2,
							const float * const __restrict__ minz1, const float * const __restrict__ minz2,
							const float * const __restrict__ maxx1, const float * const __restrict__ maxx2,
							const float * const __restrict__ maxy1, const float * const __restrict__ maxy2,
							const float * const __restrict__ maxz1, const float * const __restrict__ maxz2,
							int * const __restrict__ out
							)
		{
			alignas(32)
			int result[16]={
				// 0v0 0v1 0v2 0v3
				// 1v0 1v1 1v2 1v3
				// 2v0 2v1 2v2 2v3
				// 3v0 3v1 3v2 3v3
				0, 0, 0, 0,
				0, 0, 0, 0,
				0, 0, 0, 0,
				0, 0, 0, 0
			};

			alignas(32)
			int tileId1[16]={
					// 0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3
					partId1[0],partId1[1],partId1[2],partId1[3],
					partId1[0],partId1[1],partId1[2],partId1[3],
					partId1[0],partId1[1],partId1[2],partId1[3],
					partId1[0],partId1[1],partId1[2],partId1[3]
			};

			alignas(32)
			int tileId2[16]={
					// 0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3
					partId2[0],partId2[0],partId2[0],partId2[0],
					partId2[1],partId2[1],partId2[1],partId2[1],
					partId2[2],partId2[2],partId2[2],partId2[2],
					partId2[3],partId2[3],partId2[3],partId2[3]
			};

			alignas(32)
			float tileMinX1[16]={
					// 0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3
					minx1[0],minx1[1],minx1[2],minx1[3],
					minx1[0],minx1[1],minx1[2],minx1[3],
					minx1[0],minx1[1],minx1[2],minx1[3],
					minx1[0],minx1[1],minx1[2],minx1[3]
			};

			alignas(32)
			float tileMinX2[16]={
					// 0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3
					minx2[0],minx2[0],minx2[0],minx2[0],
					minx2[1],minx2[1],minx2[1],minx2[1],
					minx2[2],minx2[2],minx2[2],minx2[2],
					minx2[3],minx2[3],minx2[3],minx2[3]
			};

			alignas(32)
			float tileMinY1[16]={
					// 0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3
					miny1[0],miny1[1],miny1[2],miny1[3],
					miny1[0],miny1[1],miny1[2],miny1[3],
					miny1[0],miny1[1],miny1[2],miny1[3],
					miny1[0],miny1[1],miny1[2],miny1[3]
			};

			alignas(32)
			float tileMinY2[16]={
					// 0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3
					miny2[0],miny2[0],miny2[0],miny2[0],
					miny2[1],miny2[1],miny2[1],miny2[1],
					miny2[2],miny2[2],miny2[2],miny2[2],
					miny2[3],miny2[3],miny2[3],miny2[3]
			};

			alignas(32)
			float tileMinZ1[16]={
					// 0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3
					minz1[0],minz1[1],minz1[2],minz1[3],
					minz1[0],minz1[1],minz1[2],minz1[3],
					minz1[0],minz1[1],minz1[2],minz1[3],
					minz1[0],minz1[1],minz1[2],minz1[3]
			};

			alignas(32)
			float tileMinZ2[16]={
					// 0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3
					minz2[0],minz2[0],minz2[0],minz2[0],
					minz2[1],minz2[1],minz2[1],minz2[1],
					minz2[2],minz2[2],minz2[2],minz2[2],
					minz2[3],minz2[3],minz2[3],minz2[3]
			};











			alignas(32)
			float tileMaxX1[16]={
					// 0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3
					maxx1[0],maxx1[1],maxx1[2],maxx1[3],
					maxx1[0],maxx1[1],maxx1[2],maxx1[3],
					maxx1[0],maxx1[1],maxx1[2],maxx1[3],
					maxx1[0],maxx1[1],maxx1[2],maxx1[3]
			};

			alignas(32)
			float tileMaxX2[16]={
					// 0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3
					maxx2[0],maxx2[0],maxx2[0],maxx2[0],
					maxx2[1],maxx2[1],maxx2[1],maxx2[1],
					maxx2[2],maxx2[2],maxx2[2],maxx2[2],
					maxx2[3],maxx2[3],maxx2[3],maxx2[3]
			};

			alignas(32)
			float tileMaxY1[16]={
					// 0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3
					maxy1[0],maxy1[1],maxy1[2],maxy1[3],
					maxy1[0],maxy1[1],maxy1[2],maxy1[3],
					maxy1[0],maxy1[1],maxy1[2],maxy1[3],
					maxy1[0],maxy1[1],maxy1[2],maxy1[3]
			};

			alignas(32)
			float tileMaxY2[16]={
					// 0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3
					maxy2[0],maxy2[0],maxy2[0],maxy2[0],
					maxy2[1],maxy2[1],maxy2[1],maxy2[1],
					maxy2[2],maxy2[2],maxy2[2],maxy2[2],
					maxy2[3],maxy2[3],maxy2[3],maxy2[3]
			};

			alignas(32)
			float tileMaxZ1[16]={
					// 0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3
					maxz1[0],maxz1[1],maxz1[2],maxz1[3],
					maxz1[0],maxz1[1],maxz1[2],maxz1[3],
					maxz1[0],maxz1[1],maxz1[2],maxz1[3],
					maxz1[0],maxz1[1],maxz1[2],maxz1[3]
			};

			alignas(32)
			float tileMaxZ2[16]={
					// 0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3
					maxz2[0],maxz2[0],maxz2[0],maxz2[0],
					maxz2[1],maxz2[1],maxz2[1],maxz2[1],
					maxz2[2],maxz2[2],maxz2[2],maxz2[2],
					maxz2[3],maxz2[3],maxz2[3],maxz2[3]
			};




			for(int i=0;i<16;i++)
				result[i] =  (tileId1[i] < tileId2[i]);

			for(int i=0;i<16;i++)
				result[i] = result[i] &&
				intersectDim(tileMinX1[i], tileMaxX1[i], tileMinX2[i], tileMaxX2[i]) &&
				intersectDim(tileMinY1[i], tileMaxY1[i], tileMinY2[i], tileMaxY2[i]) &&
				intersectDim(tileMinZ1[i], tileMaxZ1[i], tileMinZ2[i], tileMaxZ2[i]);

			for(int i=0;i<16;i++)
				out[i]=result[i];
		};
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
			fields->mem.index.set(indexChildNodeStart,3);


			// set AABB of current (root) node
			// X

			const int indexBoundMinXFloat = fields->mem.nodeMinX.allocate(1);
			fields->mem.nodeMinX.set(indexBoundMinXFloat,fields->minCornerX);


			// Y

			const int indexBoundMinYFloat = fields->mem.nodeMinY.allocate(1);
			fields->mem.nodeMinY.set(indexBoundMinYFloat,fields->minCornerY);


			// Z

			const int indexBoundMinZFloat = fields->mem.nodeMinZ.allocate(1);
			fields->mem.nodeMinZ.set(indexBoundMinZFloat,fields->minCornerZ);


			// cell inverse width

			const int indexWidthFloat = fields->mem.nodeInvWidth.allocate(1);
			fields->mem.nodeInvWidth.set(indexWidthFloat,fields->cellWidthInv);


			// cell inverse height

			const int indexHeightFloat = fields->mem.nodeInvHeight.allocate(1);
			fields->mem.nodeInvHeight.set(indexHeightFloat,fields->cellHeightInv);


			// cell inverse depth

			const int indexDepthFloat = fields->mem.nodeInvDepth.allocate(1);
			fields->mem.nodeInvDepth.set(indexDepthFloat,fields->cellDepthInv);


			fields->mem.childNodeCount.set(fields->mem.childNodeCount.allocate(1),0);
			fields->mem.nodeCollisionMask.set(fields->mem.nodeCollisionMask.allocate(1),0);


		}


		template<typename Derived>
		inline void addParticles(const int numParticlesToAdd, Derived * const __restrict__ particles)
		{
			const int pId = fields->mem.indexParticle.allocate(numParticlesToAdd);
			const int oId = fields->mem.orderParticle.allocate(numParticlesToAdd);

			const int maxXId = fields->mem.maxX.allocate(numParticlesToAdd);
			const int maxYId = fields->mem.maxY.allocate(numParticlesToAdd);
			const int maxZId = fields->mem.maxZ.allocate(numParticlesToAdd);
			const int minXId = fields->mem.minX.allocate(numParticlesToAdd);
			const int minYId = fields->mem.minY.allocate(numParticlesToAdd);
			const int minZId = fields->mem.minZ.allocate(numParticlesToAdd);
			fields->mem.index.set(1,fields->mem.index.get(1)+numParticlesToAdd);

			for(int i=0;i<numParticlesToAdd;i++)
			{

				const IParticle<float> * const curPtr = static_cast<const IParticle<float>* const>(particles+i);


				fields->mem.indexParticle.set(pId+i,curPtr->getId());
				fields->mem.orderParticle.set(oId+i,oId+i);

				fields->mem.maxX.set(maxXId+i,curPtr->getMaxX());
				fields->mem.maxY.set(maxYId+i,curPtr->getMaxY());
				fields->mem.maxZ.set(maxZId+i,curPtr->getMaxZ());
				fields->mem.minX.set(minXId+i,curPtr->getMinX());
				fields->mem.minY.set(minYId+i,curPtr->getMinY());
				fields->mem.minZ.set(minZId+i,curPtr->getMinZ());



			}
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
				const int npdiv3 = task.nodePointer/3;
				const int numChildNodes = fields->mem.childNodeCount.get(npdiv3);

				// if this is not a leaf node, traverse all child nodes (they are sparse, so may be less than 8(8bit mask) or 64(64 bit mask))
				if(pointer<0)
				{
					// get current node's information
					const float minCornerX = fields->mem.nodeMinX.get(npdiv3);
					const float minCornerY = fields->mem.nodeMinY.get(npdiv3);
					const float minCornerZ = fields->mem.nodeMinZ.get(npdiv3);

					const float cellWidthInv = fields->mem.nodeInvWidth.get(npdiv3);
					const float cellHeightInv = fields->mem.nodeInvHeight.get(npdiv3);
					const float cellDepthInv = fields->mem.nodeInvDepth.get(npdiv3);


					const int indexStartX = std::floor((minx - minCornerX)*cellWidthInv);
					const int indexEndX = std::floor((maxx - minCornerX)*cellWidthInv);

					const int indexStartY = std::floor((miny - minCornerY)*cellHeightInv);
					const int indexEndY = std::floor((maxy - minCornerY)*cellHeightInv);

					const int indexStartZ = std::floor((minz - minCornerZ)*cellDepthInv);
					const int indexEndZ = std::floor((maxz - minCornerZ)*cellDepthInv);


					// prepare cell indicator mask (1 bit = has object, 0 bit = empty))
					uint64_t maskCellsFilled=0;
					for(int k=indexStartZ; k<=indexEndZ; k++)
					{
						if(k<0 || k>=4)
							continue;
						for(int j=indexStartY; j<=indexEndY; j++)
						{
							if(j<0 || j>=4)
								continue;
							for(int i=indexStartX; i<=indexEndX; i++)
							{
								if(i<0 || i>=4)
									continue;

								storeBit(maskCellsFilled,1,i+j*4+k*16);

							}
						}
					}


					const int nodeOffset = -pointer-1;
					for(int i=0;i<numChildNodes;i++)
					{
						// if there is possible collision (accelerated by bit mask for collisions)
						uint64_t cellMask = fields->mem.nodeCollisionMask.get((nodeOffset+i*3)/3);
						if(maskCellsFilled & cellMask)
						{
							nodesToCompute.push(NodeTask(nodeOffset+i*3));
						}
					}
				}
				else
				{
					// this is leaf node

					const int ptr = fields->mem.index.get(task.nodePointer);
					const int n = fields->mem.index.get(task.nodePointer+1);




					for(int i=0;i<n;i++)
					{
						const int index = ptr+i;

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

									fastSet.iterateSet(partId);
								}
							}
						}

				}
			}
			const int fsz = fastSet.size();
			for(int i=0;i<fsz;i++)
			{
				result.push_back(fastSet.get(i));
			}
			return result;
		}






		// visits all leaf nodes and computes nxn collision pairs
		std::vector<std::pair<int,int>> findCollisionsAll()
		{
			//std::unordered_map<int,FastUnique<int32_t, testParticleLimit>> fastSet;
			const int resetN = fields->mem.indexParticle.size();

			fields->mem.allPairsCollmapping.allocate(resetN);
			for(int i=0;i<resetN;i++)
			{
				fields->mem.allPairsCollmapping.getRef(i).reset();
			}
			fields->mem.allPairsCollmapping.reset();


			//fastSet.reserve(fields->mem.indexParticle.size());
			fields->mem.allPairsColl.reset();
			std::vector<std::pair<int,int>> result;
			//result.reserve(100000);


			const int numLeaf = fields->mem.leafOffset.size();

			for(int leaf=0;leaf<numLeaf;leaf++)
			{

				{

					const int leafOfs = fields->mem.leafOffset.get(leaf);
					const int ptr = fields->mem.index.get(leafOfs);
					const int n = fields->mem.index.get(leafOfs+1);

					alignas(32)
					int index[testParticleLimit];

					alignas(32)
					int orderId[testParticleLimit];

					alignas(32)
					int partId[testParticleLimit];

					alignas(32)
					float minx[testParticleLimit];

					alignas(32)
					float miny[testParticleLimit];

					alignas(32)
					float minz[testParticleLimit];

					alignas(32)
					float maxx[testParticleLimit];

					alignas(32)
					float maxy[testParticleLimit];

					alignas(32)
					float maxz[testParticleLimit];
					constexpr int simd = 4;
					constexpr int simd1 = simd-1;
					const int n8 = n-(n&simd1);
					for(int i=0;i<n8;i+=simd)
					{
						for(int j=0;j<simd;j++)
							index[i+j]   = ptr + i + j;
						for(int j=0;j<simd;j++)
							orderId[i+j] = fields->mem.orderParticle.get(index[i+j]);
						for(int j=0;j<simd;j++)
							partId[i+j]  = fields->mem.indexParticle.get(orderId[i+j]);
						for(int j=0;j<simd;j++)
							minx[i+j]    = fields->mem.minX.get(orderId[i+j]);
						for(int j=0;j<simd;j++)
							miny[i+j]    = fields->mem.minY.get(orderId[i+j]);
						for(int j=0;j<simd;j++)
							minz[i+j]    = fields->mem.minZ.get(orderId[i+j]);
						for(int j=0;j<simd;j++)
							maxx[i+j]    = fields->mem.maxX.get(orderId[i+j]);
						for(int j=0;j<simd;j++)
							maxy[i+j]    = fields->mem.maxY.get(orderId[i+j]);
						for(int j=0;j<simd;j++)
							maxz[i+j]    = fields->mem.maxZ.get(orderId[i+j]);

					}

					for(int i=n8;i<n;i++)
					{
							index[i]   = ptr + i;
							orderId[i] = fields->mem.orderParticle.get(index[i]);
							partId[i]  = fields->mem.indexParticle.get(orderId[i]);
							minx[i]    = fields->mem.minX.get(orderId[i]);
							miny[i]    = fields->mem.minY.get(orderId[i]);
							minz[i]    = fields->mem.minZ.get(orderId[i]);
							maxx[i]    = fields->mem.maxX.get(orderId[i]);
							maxy[i]    = fields->mem.maxY.get(orderId[i]);
							maxz[i]    = fields->mem.maxZ.get(orderId[i]);
					}

					for(int i=n;i<testParticleLimit;i++)
					{
						index[i]   = -1;
						orderId[i]   = -1;
						partId[i]   = -1;
						minx[i] = 1000000000000000000.0f;
						miny[i] = 1000000000000000000.0f;
						minz[i] = 1000000000000000000.0f;
						maxx[i] = 1000000000000000000.0f;
						maxy[i] = 1000000000000000000.0f;
						maxz[i] = 1000000000000000000.0f;
					}



					// SIMD computation (tiled computing)

					if(true)
					{
						alignas(32)
						int out[16]={
								0,0,0,0,
								0,0,0,0,
								0,0,0,0,
								0,0,0,0
						};

						for(int i=0;i<testParticleLimit-simd;i+=simd)
						{
							if(i>=n)
								break;

							FastUnique<int32_t, testParticleLimit> * map[simd] = {
									partId[i]>=0?fields->mem.allPairsCollmapping.getPtr(partId[i]):nullptr,
									partId[i+1]>=0?fields->mem.allPairsCollmapping.getPtr(partId[i+1]):nullptr,
									partId[i+2]>=0?fields->mem.allPairsCollmapping.getPtr(partId[i+2]):nullptr,
									partId[i+3]>=0?fields->mem.allPairsCollmapping.getPtr(partId[i+3]):nullptr
							};

							for(int j=i;j<testParticleLimit;j+=simd)
							{
								if(j>=n)
									break;
								// 0v0, 0v1, 0v2, 0v3,
								// 1v0, 1v1, 1v2, 1v3,
								// 2v0, 2v1, 2v2, 2v3,
								// 3v0, 3v1, 3v2, 3v3,

								comp4vs4(	partId+i, partId+j,
											minx+i, minx+j,
											miny+i, miny+j,
											minz+i, minz+j,
											maxx+i, maxx+j,
											maxy+i, maxy+j,
											maxz+i, maxz+j,
											out
								);

								for(int k=0;k<16;k++)
								{
									const int k3 = k&3;
									const int id1 = i+k3;
									const int id2 = j+(k>>2);
									if(out[k]/* && partId[id1]>=0 && partId[id2]>=0*/)
									{
										//FastUnique<int32_t, testParticleLimit>& map = fields->mem.allPairsCollmapping.getRef(partId[id1]);
										if(map[k3])
											map[k3]->insert(partId[id2]);
									}
								}
							}
						}
					}


				}

			}




			for(int i=0;i<resetN;i++)
			{
				FastUnique<int32_t, testParticleLimit>& map = fields->mem.allPairsCollmapping.getRef(i);
				const int ms = map.size();
				const int allocIdx = fields->mem.allPairsColl.allocate(ms);

				for(int j=0;j<ms;j++)
				{
					fields->mem.allPairsColl.set(allocIdx+j,std::pair<int,int>(fields->mem.indexParticle.get(i),fields->mem.indexParticle.get(map.get(j))));
				}
			}


			result.resize(fields->mem.allPairsColl.size());
			fields->mem.allPairsColl.writeTo(result);
			return result;
		}

		void buildTree()
		{



			int particleStart = fields->mem.index.get(0);
			int numParticle = fields->mem.index.get(1);

			int nodeOffset = 0;



			float minCornerX = fields->mem.nodeMinX.get(0);
			float minCornerY = fields->mem.nodeMinY.get(0);
			float minCornerZ = fields->mem.nodeMinZ.get(0);
			float cellWidthInv =  fields->mem.nodeInvWidth.get(0);
			float cellHeightInv =  fields->mem.nodeInvHeight.get(0);
			float cellDepthInv =  fields->mem.nodeInvDepth.get(0);
			float cellWidth =  1.0f/cellWidthInv;
			float cellHeight =  1.0f/cellHeightInv;
			float cellDepth =  1.0f/cellDepthInv;
			int ctr=0;

			int maxNodeOffset = 3;

			while(nodeOffset <= maxNodeOffset)
			{
				ctr++;

				int ctrTmp[64]={0,0,0,0,0,0,0,0,
						0,0,0,0,0,0,0,0,
						0,0,0,0,0,0,0,0,
						0,0,0,0,0,0,0,0,
						0,0,0,0,0,0,0,0,
						0,0,0,0,0,0,0,0,
						0,0,0,0,0,0,0,0,
						0,0,0,0,0,0,0,0};


				// if child node pointer not set up


				if(fields->mem.index.get(nodeOffset+2)<nodeOffset && fields->mem.index.get(fields->mem.index.get(nodeOffset+2)+2)>=0)
				{
					fields->mem.index.set(fields->mem.index.get(nodeOffset+2)+2,-(nodeOffset+1));
				}


				int childNodeCount = 0;

				if(numParticle > testParticleLimit)
				{


					{

						for(int zz = 0; zz<4; zz++)
							for(int yy = 0; yy<4; yy++)
								for(int xx = 0; xx<4; xx++)
								{
									// allocate node
									const int index0 = xx+yy*4+zz*16;

									fields->mem.orderTmp[index0].reset();
									fields->mem.orderTmp[index0].allocate(numParticle);
								}
					}


					{




						for(int ii=0;ii<numParticle;ii++)
						{

							const int orderParticle = fields->mem.orderParticle.get(particleStart+ii);
							const float& minx = fields->mem.minX.getRef(orderParticle);
							const float& miny = fields->mem.minY.getRef(orderParticle);
							const float& minz = fields->mem.minZ.getRef(orderParticle);
							const float& maxx = fields->mem.maxX.getRef(orderParticle);
							const float& maxy = fields->mem.maxY.getRef(orderParticle);
							const float& maxz = fields->mem.maxZ.getRef(orderParticle);

							const int indexStartX = std::floor((minx - minCornerX)*cellWidthInv);
							const int indexEndX = std::floor((maxx - minCornerX)*cellWidthInv);

							const int indexStartY = std::floor((miny - minCornerY)*cellHeightInv);
							const int indexEndY = std::floor((maxy - minCornerY)*cellHeightInv);

							const int indexStartZ = std::floor((minz -minCornerZ)*cellDepthInv);
							const int indexEndZ = std::floor((maxz - minCornerZ)*cellDepthInv);

							// prepare cell indicator mask (1 bit = has object, 0 bit = empty))

							for(int k=indexStartZ; k<=indexEndZ; k++)
							{
								if(k<0 || k>=4)
									continue;
								for(int j=indexStartY; j<=indexEndY; j++)
								{
									if(j<0 || j>=4)
										continue;
									for(int i=indexStartX; i<=indexEndX; i++)
									{
										if(i<0 || i>=4)
											continue;


										const int index0 = i+j*4+k*16;

										fields->mem.orderTmp[index0].set(ctrTmp[index0],orderParticle);
										ctrTmp[index0]++;
									}
								}
							}
						}


					}


					// add all particles in order (from first child node to last child node)
					childNodeCount=0;
					for(int zz = 0; zz<4; zz++)
						for(int yy = 0; yy<4; yy++)
							for(int xx = 0; xx<4; xx++)
							{
								const int index0 = xx+yy*4+zz*16;
								const int sz = ctrTmp[index0];


								if(sz>0)
								{
									childNodeCount++;

									const int nodeIndexOfs         = fields->mem.index.allocate(3);
									const int particleStartCur     = nodeIndexOfs;
									const int numParticleCur       = nodeIndexOfs+1;
									const int childNodeStartCur    = nodeIndexOfs+2;

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

									fields->mem.nodeInvWidth.set(nodeInvWidthFloat,cellWidthInv*4.0f);
									fields->mem.nodeInvHeight.set(nodeInvHeightFloat,cellHeightInv*4.0f);
									fields->mem.nodeInvDepth.set(nodeInvDepthFloat,cellDepthInv*4.0f);


									const int nodeMaskIndex = fields->mem.nodeCollisionMask.allocate(1);
									uint64_t nodeMask = 0;
									storeBit(nodeMask,1,index0);
									fields->mem.nodeCollisionMask.set(nodeMaskIndex,nodeMask);

									//const int allocOffset = fields->mem.indexParticle.allocate(sz);
									const int allocOffset = fields->mem.orderParticle.allocate(sz);


									//fields->mem.indexParticle.readFrom(fields->mem.idTmp[index0],0,allocOffset,sz);
									fields->mem.orderParticle.readFrom(fields->mem.orderTmp[index0],0,allocOffset,sz);


									fields->mem.index.set(particleStartCur,allocOffset);
									fields->mem.index.set(numParticleCur,sz);
									fields->mem.index.set(childNodeStartCur,nodeOffset);


									maxNodeOffset=particleStartCur;
								}



							}

					fields->mem.childNodeCount.set(nodeOffset/3,childNodeCount);

				}
				else
				{
					fields->mem.childNodeCount.set(nodeOffset/3,0);
					const int idx = fields->mem.leafOffset.allocate(1);
					fields->mem.leafOffset.set(idx,nodeOffset);
				}

				nodeOffset += 3;
				numParticle=0;
				if(nodeOffset <= maxNodeOffset)
				{
					particleStart = fields->mem.index.get(nodeOffset);
					numParticle = fields->mem.index.get(nodeOffset+1);


					minCornerX = fields->mem.nodeMinX.get(nodeOffset/3);
					minCornerY = fields->mem.nodeMinY.get(nodeOffset/3);
					minCornerZ = fields->mem.nodeMinZ.get(nodeOffset/3);
					cellWidthInv =  fields->mem.nodeInvWidth.get(nodeOffset/3);
					cellHeightInv =  fields->mem.nodeInvHeight.get(nodeOffset/3);
					cellDepthInv =  fields->mem.nodeInvDepth.get(nodeOffset/3);
					cellWidth =  1.0f/cellWidthInv;
					cellHeight =  1.0f/cellHeightInv;
					cellDepth =  1.0f/cellDepthInv;
				}
			}




		}

	private:
		std::shared_ptr<AdaptiveGridV2Fields> fields;
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
