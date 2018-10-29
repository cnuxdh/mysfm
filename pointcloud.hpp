

#ifndef POINT_CLOUD_HPP
#define POINT_CLOUD_HPP


struct RegisterPtCloudAndImageError
{
	RegisterPtCloudAndImageError(double observed_x, double observed_y)
		: observed_x(observed_x), observed_y(observed_y) {}

	template <typename T>
	bool operator()(const T* const params, 
		T* residuals) const   
	{

	}

	double observed_x, observed_y;

};
















#endif
