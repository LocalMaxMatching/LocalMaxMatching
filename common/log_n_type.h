
namespace dipl
{

template <typename ElementType>
struct LogNType
{
	typedef int Type;
};

template <>
struct LogNType<int>
{
	typedef char Type;
};

template <>
struct LogNType<unsigned int>
{
	typedef unsigned char Type;
};

template <>
struct LogNType<long>
{
	typedef char Type;
};

template <>
struct LogNType<unsigned long>
{
	typedef unsigned char Type;
};

} // end namespace dipl
