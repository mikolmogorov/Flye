//(c) 2019 by Authors
//This file is a part of the Flye program.
//Released under the BSD license (see LICENSE file)

//Big Functional Container - a substitute for vector
//for very large arrays (up to several billion elements).
//Basically, it is a vector of vectors - so we have 
//a list of (relatively big) chunks instead of one
//gigantic memory chunk (several Gb) - allocating which might be
//difficult for OS, cause memory fragmentation and slowdown.
//Includes memory pool of chunks that could be shared by multuple threads

#pragma once

#include <vector>
#include <mutex>

template <class T, int ChunkSize = 1024 * 1024>
class ChunkPool
{
public:
	ChunkPool(){}
	ChunkPool(const ChunkPool& other) = delete;
	ChunkPool(ChunkPool&& other) = delete;
	ChunkPool& operator=(ChunkPool& other) = delete;

	~ChunkPool()

	{
		for (T* chunk : _freeChunks) delete[] chunk;
		_freeChunks.clear();
	}

	T* getChunk()
	{
		std::lock_guard<std::mutex> lock(_chunkMutex);
		if (!_freeChunks.empty())
		{
			T* ref = _freeChunks.back();
			_freeChunks.pop_back();
			return ref;
		}
		else
		{
			return new T[ChunkSize];
		}
	}

	void returnChunk(T* chunk)
	{
		std::lock_guard<std::mutex> lock(_chunkMutex);
		_freeChunks.push_back(chunk);
	}

	size_t numberChunks() {return _freeChunks.size();}
	
private:
	std::mutex 		_chunkMutex;
	std::vector<T*> _freeChunks;
};


template <class T, size_t ChunkSize = 1024 * 1024>
class BFContainer
{
public:
	BFContainer(ChunkPool<T, ChunkSize>& chunkPool): 
		_pool(chunkPool),
		_size(0), _lastChunkOffset(0) 
	{
		_chunks.push_back(_pool.getChunk());
	}

	BFContainer(const BFContainer& other) = delete;
	BFContainer(BFContainer&& other) = delete;
	BFContainer& operator=(const BFContainer& other) = delete;

	~BFContainer()
	{
		for (T* chunk : _chunks)
		{
			_pool.returnChunk(chunk);
		}
		_chunks.clear();
	}

	void push_back(T&& elem)
	{
		_chunks.back()[_lastChunkOffset] = elem;
		++_lastChunkOffset;
		++_size;
		if (_lastChunkOffset == ChunkSize)
		{
			_chunks.push_back(_pool.getChunk());
			_lastChunkOffset = 0;
		}
	}

	template <typename ... Args>
	void emplace_back(Args&& ... args)
	{
		new (_chunks.back() + _lastChunkOffset) T(std::forward<Args>(args)...);
		++_lastChunkOffset;
		++_size;
		if (_lastChunkOffset == ChunkSize)
		{
			_chunks.push_back(_pool.getChunk());
			_lastChunkOffset = 0;
		}
	}

	size_t size() {return _size;}

	T& operator[](size_t index)
	{
		const size_t chunkId = index / ChunkSize;
		const size_t chunkOffset = index - chunkId * ChunkSize;
		return _chunks[chunkId][chunkOffset];
	}

	void clear()
	{
	}

	class BFIterator: public std::iterator<std::random_access_iterator_tag, T>
	{
	private:
		T** _map;
		T*  _chunkStart;
		T*  _chunkCur;
		T*  _chunkEnd;

    public:
		typedef ptrdiff_t difference_type;

        BFIterator(T** map=nullptr, T* chunkStart=nullptr,
				   T* chunkCur=nullptr, T* chunkEnd=nullptr): 
			_map(map), _chunkStart(chunkStart), 
			_chunkCur(chunkCur), _chunkEnd(chunkEnd) {}

        BFIterator& operator=(const BFIterator &rhs) 
		{
			_map = rhs._map; 
			_chunkStart = rhs._chunkStart; 
			_chunkCur = rhs._chunkCur; 
			_chunkEnd = rhs._chunkEnd; 
			return *this;
		}

		void setNode(T** node)
		{
			_map = node;
			_chunkStart = *_map;
			_chunkEnd = _chunkStart + ChunkSize;
		}

        BFIterator& operator+=(difference_type n) 
		{
			const difference_type offset = n + (_chunkCur - _chunkStart);
			if (offset >= 0 && offset < difference_type(ChunkSize))
			{
				_chunkCur += n;
			}
			else
			{
				const difference_type nodeOffset = 
					(offset > 0) ? offset / difference_type(ChunkSize) :
			   					   -difference_type((-offset - 1) / ChunkSize) - 1;
				setNode(_map + nodeOffset);
				_chunkCur = _chunkStart + (offset - nodeOffset * ChunkSize);
			}
			return *this;
		}

        BFIterator& operator-=(difference_type n) 
		{
			return *this += -n;
		}

		BFIterator operator+(difference_type n) const
		{
			auto tmp = *this;
			return tmp += n;
		}

		friend BFIterator operator+(difference_type lhs, const BFIterator& rhs) 
		{
			return rhs + lhs;
		}
		
        BFIterator operator-(difference_type n) const
		{
			auto tmp = *this;
			return tmp -= n;
		}

        difference_type operator-(const BFIterator& rhs) const
		{
			return difference_type(ChunkSize) * (_map - rhs._map - 1) + 
						(_chunkCur - _chunkStart) +
						(rhs._chunkEnd - rhs._chunkCur);
		}

		BFIterator& operator++() 
		{
			++_chunkCur;
			if (_chunkCur == _chunkEnd)
	  		{
	    		setNode(_map + 1);
	    		_chunkCur = _chunkStart;
	  		}
			return *this;
		}

        BFIterator operator++(int)
		{
			auto tmp(*this); 
			++(*this);
			return tmp;
		}

        BFIterator& operator--() 
		{
			if (_chunkCur == _chunkStart)
	  		{
	    		setNode(_map - 1);
	    		_chunkCur = _chunkEnd;
	  		}
			--_chunkCur;
			return *this;
		}

        BFIterator operator--(int)
		{
			auto tmp(*this); 
			--(*this);
			return tmp;
		}

		T& operator[](difference_type n) {return *(*this += n);}

		T& operator*() const {return *_chunkCur;}
        T* operator->() const {return _chunkCur;}

    	bool operator==(const BFIterator& rhs) const
		{
			return _chunkCur == rhs._chunkCur;
		}
        bool operator!=(const BFIterator& rhs) const
		{
			return !(*this == rhs);
		}
        bool operator< (const BFIterator& rhs) const
		{
			return _map == rhs._map ? _chunkCur < rhs._chunkCur : _map < rhs._map;
		}
        bool operator> (const BFIterator& rhs) const
		{
			return rhs < *this;
		}
        bool operator<=(const BFIterator& rhs) const
		{
			return _map == rhs._map ? _chunkCur <= rhs._chunkCur : _map < rhs._map;
		}
        bool operator>=(const BFIterator& rhs) const
		{
			return rhs <= *this;
		}
	};

	BFIterator begin() {return BFIterator(&_chunks.front(), _chunks.front(), 
										  _chunks.front(), 
										  _chunks.front() + ChunkSize);}

	BFIterator end() {return BFIterator(&_chunks.back(), _chunks.back(), 
										_chunks.back() + _lastChunkOffset, 
										_chunks.back() + ChunkSize);}

private:

	ChunkPool<T, ChunkSize>& 	_pool;
	std::vector<T*> 			_chunks;
	size_t 						_size;
	size_t 						_lastChunkOffset;
};
