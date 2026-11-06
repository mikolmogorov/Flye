//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include <cstdlib>
#include <stdexcept>
#include <cassert>
#include <utility>
#include <algorithm>

template <class T>
class Matrix
{
public:
	Matrix(): _rows(0), _cols(0), _capacity(0), _data(nullptr) {}
	Matrix(const Matrix& other):
		Matrix(other._rows, other._cols)
	{
		std::copy(other._data, other._data + _rows * _cols, _data);
	}
	Matrix(Matrix&& other) noexcept:
		Matrix()
	{
		std::swap(_cols, other._cols);
		std::swap(_rows, other._rows);
		std::swap(_capacity, other._capacity);
		std::swap(_data, other._data);
	}
	Matrix& operator=(Matrix && other) noexcept
	{
		std::swap(_cols, other._cols);
		std::swap(_rows, other._rows);
		std::swap(_capacity, other._capacity);
		std::swap(_data, other._data);
		return *this;
	}
	Matrix& operator=(const Matrix& other)
	{
		if (this != &other)
		{
			this->resize(other._rows, other._cols);
			std::copy(other._data, other._data + _rows * _cols, _data);
		}
		return *this;
	}

	Matrix(size_t rows, size_t cols, T val = 0):
		_rows(rows), _cols(cols), _capacity(rows * cols)
	{
		if (!rows || !cols)
			throw std::runtime_error("Zero matrix dimension");
		_data = new T[_capacity];
		for (size_t i = 0; i < _capacity; ++i) _data[i] = val;
	}
	void resize(size_t rows, size_t cols)
	{
		if (!rows || !cols)
			throw std::runtime_error("Zero matrix dimension");
		size_t needed = rows * cols;
		if (needed > _capacity)
		{
			delete[] _data;
			_capacity = needed;
			_data = new T[_capacity];
		}
		_rows = rows;
		_cols = cols;
		_data[0] = T(0);
	}
	~Matrix()
	{
		if (_data) delete[] _data;
	}

	T* rowPtr(size_t row) {return _data + row * _cols;}
	const T* rowPtr(size_t row) const {return _data + row * _cols;}
	T& at(size_t row, size_t col) {return _data[row * _cols + col];}
	const T& at(size_t row, size_t col) const {return _data[row * _cols + col];}
	size_t nrows() const {return _rows;}
	size_t ncols() const {return _cols;}

private:
	size_t _rows;
	size_t _cols;
	size_t _capacity;
	T* _data;
};
