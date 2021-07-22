#pragma once
#include <algorithm>
#include <sstream>
#include <vector>
#include <fstream>
#include <map>

template<class T> class array3d {
	std::unique_ptr<T*> array = nullptr;
	std::size_t x_, y_, z_;
	std::size_t size_;
public:
	array3d<T>(int x1 = 0, int y1 = 0, int z1 = 0);
	array3d<T>(int x1, int y1, int z1, T val);
	array3d<T>(const array3d<T>&);
	array3d<T>& operator=(const array3d<T>&);
	~array3d();
	T& operator() (int i, int j, int k) { return (*array)[i + j * x_ + k * x_ * y_]; };
	T   operator() (int i, int j, int k)const { return (*array)[i + j * x_ + k * x_ * y_]; };
	T& operator() (int i) { return (*array)[i]; }
	T operator() (int i) const { return (*array)[i]; }
	T* mass(void) { return *array; }
	const T* mass(void) const { return *array; }
	void resize(int x_, int y_, int z_);
	std::size_t x()const { return x_; }
	std::size_t y()const { return y_; }
	std::size_t z()const { return z_; }
	std::size_t size()const { return size_; }
	T max()
	{
		T max = (*array)[0];

		for (int i = 1; i < size_; ++i)
		{
			if ((*array)[i] > max) max = (*array)[i];
		}
		return max;
	}
};

template<class T> std::istream& operator>>(std::istream& is, array3d<T>& a);
template<class T> std::ostream& operator>>(std::ostream& is, const array3d<T>& a);

//implementation

template<class T> std::istream& operator>>(std::istream& is, array3d<T>& a) {
	int x, y, z;
	is >> x >> y >> z;
	a.resize(x, y, z);
	for (int k = 0; k < a.z(); ++k)
		for (int j = 0; j < a.y(); ++j)
			for (int i = 0; i < a.x(); ++i)
				is >> a(i, j, k);
	return is;
}

template<class T> std::ostream& operator<<(std::ostream& os, const array3d<T>& a) {
	os << a.x() << ' ' << a.y() << ' ' << a.z() << '\n';
	for (int k = 0; k < a.z(); ++k) {
		for (int j = 0; j < a.y(); ++j) {
			for (int i = 0; i < a.x(); ++i)
				os << a(i, j, k) << ' ';
			os << '\n';
		}
		os << "\n\n";
	}
	return os;
}

template<class T> array3d<T>::array3d(int x1, int y1, int z1) {
	x_ = x1; y_ = y1; z_ = z1; size_ = x_ * y_ * z_;
	if (size_ != 0) {
		array.release();
		array = std::make_unique<T*>(new T[size_]);
	}
};

template<class T> array3d<T>::array3d(int x1, int y1, int z1, T val) {
	x_ = x1; y_ = y1; z_ = z1; size_ = x_ * y_ * z_;
	if (size_ != 0) {
		array.release();
		array = std::make_unique<T*>(new T[size_]);
	}
	for (size_t i = 0; i < size_; ++i) (*array)[i] = val;
};

template<class T> array3d<T>::array3d(const array3d<T>& a) {
	x_ = a.x();
	y_ = a.y();
	z_ = a.z();
	size_ = a.size();
	if (size_ != 0) {
		array.release();
		array = std::make_unique<T*>(new T[size_]);
		std::copy(*a.array, *a.array + a.size(), *array);
	}
}

template<class T> array3d<T>& array3d<T>::operator=(const array3d<T>& a) {
	resize(a.x(), a.y(), a.z());
	if (size_ != 0) std::copy(*a.array, *a.array + a.size(), *array);
	return *this;
}

template<class T> array3d<T>::~array3d<T>() {
	delete[] * array;
};

template<class T> void array3d<T>::resize(int x1, int y1, int z1) {
	x_ = x1; y_ = y1; z_ = z1; size_ = x_ * y_ * z_;
	if (array)
		delete[] * array;
	array = std::make_unique<T*>(new T[size_]);
}

array3d<int> read_am_int(const char* fname, int x, int y, int z) {
	typedef unsigned short T;
	const int S = sizeof(T);
	T* a = new T[x * y * z * S];

	FILE* f = fopen(fname, "r");
	int gg = fread(a, S, x * y * z, f);

	array3d<int> A(x, y, z);
	for (size_t k = 0; k < A.z(); ++k)
		for (size_t j = 0; j < A.y(); ++j)
			for (size_t i = 0; i < A.x(); ++i) {
				unsigned char* ch1 = (unsigned char*)(a + (i + j * y + k * x * y));
				unsigned char* ch2 = ch1 + 1;
				A(i, j, k) = 256 * (*ch1) + (*ch2);
			}

	fclose(f);
	delete[] a;
	return A;
}

array3d<int> read_am_char(const char* fname, int x, int y, int z) {
	typedef unsigned char T;
	const int S = sizeof(T);

	array3d<unsigned char> A(x, y, z);
	FILE* f = fopen(fname, "r");
	int gg = fread(A.mass(), S, x * y * z, f);
	fclose(f);

	array3d<int> B(x, y, z);
	if (B.size() != 0) std::copy_n(B.mass(), B.size(), A.mass());
	return B;
}

void array3d2paraview(const array3d<int>& A, const std::string& fname) {
	FILE* f = fopen(fname.c_str(), "w");
	fprintf(f, "NDims = 3\nDimSize = %d %d %d\n", A.x(), A.y(), A.z());
	fputs("ElementSize = 1.0 1.0 1.0\nElementSpacing = 1.0 1.0 1.0\nElementType = MET_INT\nElementByteOrderMSB = False\nElementDataFile = local\n", f);
	fwrite(A.mass(), sizeof(int), A.x() * A.y() * A.z(), f);
	fclose(f);
}

void array3d2paraview(const array3d<unsigned char>& A, const std::string& fname) {
	FILE* f = fopen(fname.c_str(), "w");
	fprintf(f, "NDims = 3\nDimSize = %d %d %d\n", A.x(), A.y(), A.z());
	fputs("ElementSize = 1.0 1.0 1.0\nElementSpacing = 1.0 1.0 1.0\nElementType = MET_UCHAR\nElementByteOrderMSB = False\nElementDataFile = local\n", f);
	fwrite(A.mass(), sizeof(unsigned char), A.x() * A.y() * A.z(), f);
	fclose(f);
}

void array3d2paraview(const array3d<unsigned short>& A, const std::string& fname) {
	FILE* f = fopen(fname.c_str(), "w");
	fprintf(f, "NDims = 3\nDimSize = %d %d %d\n", A.x(), A.y(), A.z());
	fputs("ElementSize = 1.0 1.0 1.0\nElementSpacing = 1.0 1.0 1.0\nElementType = MET_USHORT\nElementByteOrderMSB = False\nElementDataFile = local\n", f);
	fwrite(A.mass(), sizeof(unsigned short), A.x() * A.y() * A.z(), f);
	fclose(f);
}
//Don't work correct if numbers different
void array3d2paraview(const array3d<double>& A, const std::string& fname) {
	FILE* f = fopen(fname.c_str(), "w");
	fprintf(f, "NDims = 3\nDimSize = %d %d %d\n", A.x(), A.y(), A.z());
	fputs("ElementSize = 1.0 1.0 1.0\nElementSpacing = 1.0 1.0 1.0\nElementType = MET_DOUBLE\nElementByteOrderMSB = False\nElementDataFile = local\n", f);
	fwrite(A.mass(), sizeof(double), A.size(), f);
	fclose(f);
}

namespace
{
	template<class T>
	void SwapEnd(T& var)
	{
		unsigned char* varArray = reinterpret_cast<unsigned char*>(&var);
		for (size_t i = 0; i < static_cast<size_t>(sizeof(var) / 2); ++i)
			std::swap(varArray[sizeof(var) - 1 - i], varArray[i]);
	}
}

//Write filename with extension. You can save data only in int, float and double! 
template<class T>
bool array3d2vtk(const array3d<T>& A, const std::string& fname)
{
	if (fname.empty())
	{
		return false;
	}
	std::string data_type = typeid(T).name();
	if (data_type != "double" && data_type != "float" && data_type != "int" && data_type != "char" && data_type != "unsigned short" && data_type != "unsigned char")
	{
		return false;
	}
	if (data_type == "unsigned short")
	{
		data_type = "unsigned_short";
	}
	if (data_type == "unsigned char")
	{
		data_type = "unsigned_char";
	}
	std::ofstream file;
	file.open(fname, std::ios::out | std::ios::binary);
	if (!file.is_open()) { return false; }
	file << "# vtk DataFile Version 5.1" << std::endl
		<< "vtk output" << std::endl
		<< "BINARY" << std::endl
		<< "DATASET STRUCTURED_POINTS" << std::endl
		<< "DIMENSIONS " << A.x() << " " << A.y() << " " << A.z() << std::endl
		<< "SPACING 1 1 1" << std::endl
		<< "ORIGIN 0 0 0" << std::endl
		<< "POINT_DATA " << A.size() << std::endl
		<< "SCALARS MetaImage " << data_type << std::endl
		<< "LOOKUP_TABLE default" << std::endl;
	for (size_t i = 0; i < A.size(); ++i)
	{
		T tmp = A.mass()[i];
		SwapEnd(tmp);
		file.write(reinterpret_cast<char*>(&tmp), sizeof(T));
	}
	file << std::endl;
	file << "METADATA" << std::endl
		<< "INFORMATION 0" << std::endl;

	file.close();
	return true;
}

template<class T>
bool vtk2array3d(const std::string& fname, array3d<T>& data)
{
	if (fname.empty())
	{
		return false;
	}
	std::ifstream file;
	file.open(fname, std::ios::in | std::ios::binary);
	if (!file.is_open())
	{
		return false;
	}
	std::string tmp;
	for (int i = 0; i < 11; ++i)
	{
		file >> tmp;
	}
	int x, y, z;
	file >> x >> y >> z;
	for (int i = 0; i < 15; ++i)
	{
		file >> tmp;
	}
	file.read(tmp.data(), 1);
	data = array3d<T>(x, y, z);
	for (size_t i = 0; i < data.size(); ++i)
	{
		T tmp_data;
		file.read(reinterpret_cast<char*>(&tmp_data), sizeof(T));
		SwapEnd(tmp_data);
		data(i) = tmp_data;
	}
}

bool am2vtk(const std::string& source_fname, const std::string& save_fname)
{
	if (source_fname.empty() || save_fname.empty())
	{
		return false;
	}
	std::string load_name = source_fname + ".am";
	std::string save_name = save_fname + ".vtk";
	std::ofstream save_file;
	std::ifstream load_file;
	load_file.open(load_name, std::ios::out | std::ios::binary);
	save_file.open(load_name, std::ios::in | std::ios::binary);
	if (!load_file.is_open() || !save_file.is_open())
	{
		return false;
	}
	std::string tmp;
	for (int i = 0; i < 10; ++i) {
		load_file >> tmp;
	}
	return false;
}

namespace {
	void ToUpper(std::string& input)
	{
		std::for_each(std::begin(input), std::end(input), [](char& c) {
			c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
			});
	}

	template<class T>
	bool mhd2Anyarray3d(const std::string& source, array3d<T>& output)
	{
		std::ifstream load_file;
		load_file.open(source, std::ios::out | std::ios::binary);
		if (!load_file.is_open())
		{
			return false;
		}
		std::string tmp;
		int x, y, z;

		for (int i = 0; i < 5; ++i)
		{
			load_file >> tmp;
		}
		load_file >> x >> y >> z;
		for (int i = 0; i < 11; ++i)
		{
			load_file >> tmp;
		}
		if (tmp.compare("ElementType"))
		{
			return false;//Element type didn't find
		}
		load_file >> tmp >> tmp;
		std::string type_ = typeid(T).name();
		if (!type_.compare("unsigned short"))
		{
			type_ = "USHORT";
		}
		if (!type_.compare("unsigned char"))
		{
			type_ = "UCHAR";
		}
		std::string type_name = type_;
		ToUpper(type_name);
		type_name = "MET_" + type_name;
		if (tmp.compare(type_name.c_str()))
		{
			return false;//wrong type
		}
		for (int i = 0; i < 6; ++i)
			load_file >> tmp;
		load_file.read(tmp.data(), 1);
		output = array3d<T>(x, y, z);
		for (size_t i = 0; i < output.size(); ++i)
		{
			T tmp_data;
			load_file.read(reinterpret_cast<char*>(&tmp_data), sizeof(T));
			//SwapEnd(tmp_data);
			output(i) = tmp_data;
		}
		return true;
	}
}
inline bool mhd2array3d(const std::string& source, array3d<unsigned char>& output)
{
	return mhd2Anyarray3d(source, output);
}
inline bool mhd2array3d(const std::string& source, array3d<char>& output)
{
	return mhd2Anyarray3d(source, output);
}

inline bool mhd2array3d(const std::string& source, array3d<int>& output)
{
	return mhd2Anyarray3d(source, output);
}

inline bool mhd2array3d(const std::string& source, array3d<float>& output)
{
	return mhd2Anyarray3d(source, output);
}

inline bool mhd2array3d(const std::string& source, array3d<double>& output)
{
	return mhd2Anyarray3d(source, output);
}

inline bool mhd2array3d(const std::string& source, array3d<unsigned short>& output)
{
	return mhd2Anyarray3d(source, output);
}

template<class T>
void fromMhdToVtk(const std::string& source, const std::string& destination)
{
	array3d<T> temp;
	mhd2array3d(source, temp);
	array3d<T> temp2 = array3d<T>(temp.x(), temp.y(), temp.z());
	for (int i = 0; i < temp.x(); ++i) {
		for (int j = 0; j < temp.y(); ++j) {
			for (int k = 0; k < temp.z(); ++k) {
				temp2(i, j, k) = temp(i, j, k);
			}
		}
	}
	array3d2vtk(temp2, destination);
}


template<class T>
void clipVtk(const std::string& source, const std::string& destination, unsigned int size_x, unsigned int size_y, unsigned int size_z)
{
	array3d<T> temp;
	vtk2array3d(source, temp);
	array3d<T> temp2 = array3d<T>(size_x, size_y, size_z);
	for (int i = 0; i < size_x; ++i) {
		for (int j = 0; j < size_y; ++j) {
			for (int k = 0; k < size_z; ++k) {
				temp2(i, j, k) = temp(i, j, k);
			}
		}
	}
	array3d2vtk(temp2, destination);
}

template<typename T>
void clipCylinder(array3d<T>* data, double radius)
{
	if (data == nullptr)
	{
		return;
	}
	for (int i = 0; i < data->x(); ++i) {
		for (int j = 0; j < data->y(); ++j) {
			for (int k = 0; k < data->z(); ++k)
			{
				if (sqrt(pow(radius - i, 2) + pow(radius - j, 2)) > radius)
				{
					(*data)(i, j, k) = 0;
				}
			}
		}
	}
}
