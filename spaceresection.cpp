#include "spaceresection.h"

SpaceResection &SpaceResection::setInitialValue(const std::vector<double> &elems)
{
	auto i = elems.cbegin();

	m_fx = *i++;
	m_fy = *i++;
	m_x0 = *i++;
	m_y0 = *i++;

	m_k1 = *i++;
	m_k2 = *i++;
	m_p1 = *i++;
	m_p2 = *i++;

	m_width = *i++;
	m_height = *i++;

	m_X = *i++;
	m_Y = *i++;
	m_Z = *i++;
	m_phi = angleToRadian(*i++);
	m_omega = angleToRadian(*i++);
	m_kappa = angleToRadian(*i);

	return *this;
}

bool SpaceResection::importImagePoints(const std::string &file)
{
	std::fstream fin(file, std::ios::in);

	if (!fin)
		return false;

	size_t idex = 0;
	double x = 0, y = 0;
	while (!fin.eof())
	{
		fin >> idex;
		fin >> x;
		fin >> y;

		Point2D pt(idex, x, y);
		m_imagePoints.push_back(pt);
	}
	fin.close();

	return true;
}

bool SpaceResection::importControlPoints(const std::string &file)
{
	std::fstream fin(file, std::ios::in);

	if (!fin)
		return false;

	size_t idex = 0;
	double X = 0, Y = 0, Z = 0;
	while (!fin.eof())
	{
		fin >> idex;
		fin >> X;
		fin >> Y;
		fin >> Z;

		Point3D pt(idex, X, Y, Z);
		m_controlPoints.push_back(pt);			
	}
	fin.close();

	return true;
}

const std::vector<double> &SpaceResection::slove(const std::string &file)
{
	std::vector<Point3D> controlPoints;

	for (auto i = m_imagePoints.cbegin(); i != m_imagePoints.cend(); ++i)
		for (auto j = m_controlPoints.cbegin(); j != m_controlPoints.cend(); ++j)
		{
			if (i->idex == j->idex)
			{
				controlPoints.push_back(*j);
				break;
			}

			if (j == m_controlPoints.cend() - 1)
			{
				std::clog << "[error] image points and control points do not match." << std::endl;
				std::clog << "image point " << i->idex << " can not find its matched control point." << std::endl;

				system("pause");
				exit(0);
			}
		}

	m_controlPoints = controlPoints;
	controlPoints.clear();

	cv::Mat_<double> A(cv::Size(14, m_imagePoints.size() * 2), 0);
	cv::Mat_<double> L(cv::Size(1, m_imagePoints.size() * 2), 0);

	size_t iterate_times = 0;
	while (true)
	{
		double a1 = cos(m_phi)*cos(m_kappa) - sin(m_phi)*sin(m_omega)*sin(m_kappa);
		double a2 = -cos(m_phi)*sin(m_kappa) - sin(m_phi)*sin(m_omega)*cos(m_kappa);
		double a3 = -sin(m_phi)*cos(m_omega);
		double b1 = cos(m_omega)*sin(m_kappa);
		double b2 = cos(m_omega)*cos(m_kappa);
		double b3 = -sin(m_omega);
		double c1 = sin(m_phi)*cos(m_kappa) + cos(m_phi)*sin(m_omega)*sin(m_kappa);
		double c2 = -sin(m_phi)*sin(m_kappa) + cos(m_phi)*sin(m_omega)*cos(m_kappa);
		double c3 = cos(m_phi)*cos(m_omega);

		for (size_t i = 0; i < A.rows; i += 2)
		{
			double X_ = a1*(m_controlPoints[i / 2].getX() - m_X) + b1*(m_controlPoints[i / 2].getY() - m_Y) + c1*(m_controlPoints[i / 2].getZ() - m_Z);
			double Y_ = a2*(m_controlPoints[i / 2].getX() - m_X) + b2*(m_controlPoints[i / 2].getY() - m_Y) + c2*(m_controlPoints[i / 2].getZ() - m_Z);
			double Z_ = a3*(m_controlPoints[i / 2].getX() - m_X) + b3*(m_controlPoints[i / 2].getY() - m_Y) + c3*(m_controlPoints[i / 2].getZ() - m_Z);

			Point2D pt_uv = m_imagePoints[i / 2], pt_xy;
			pt_xy.x = pt_uv.x - m_x0;
			pt_xy.y = pt_uv.y - m_y0;

			double r2 = pt_xy.x*pt_xy.x + pt_xy.y*pt_xy.y;
			A.at<double>(i, 0) = 1 / Z_*(a1*m_fx + a3*pt_xy.x);
			A.at<double>(i, 1) = 1 / Z_*(b1*m_fx + b3*pt_xy.x);
			A.at<double>(i, 2) = 1 / Z_*(c1*m_fx + c3*pt_xy.x);
			A.at<double>(i, 3) = pt_xy.y*sin(m_omega) - (pt_xy.x / m_fx*(pt_xy.x*cos(m_kappa) - pt_xy.y*sin(m_kappa)) + m_fx*cos(m_kappa))*cos(m_omega);
			A.at<double>(i, 4) = -m_fx*sin(m_kappa) - pt_xy.x / m_fx*(pt_xy.x*sin(m_kappa) + pt_xy.y*cos(m_kappa));
			A.at<double>(i, 5) = pt_xy.y;
			A.at<double>(i, 6) = pt_xy.x / m_fx;
			A.at<double>(i, 7) = 0;
			A.at<double>(i, 8) = 1;
			A.at<double>(i, 9) = 0;
			A.at<double>(i, 10) = pt_xy.x*r2;
			A.at<double>(i, 11) = pt_xy.x*r2*r2;
			A.at<double>(i, 12) = r2 + 2 * pt_xy.x*pt_xy.x;
			A.at<double>(i, 13) = 2 * pt_xy.x*pt_xy.y;

			A.at<double>(i + 1, 0) = 1 / Z_*(a2*m_fy + a3*pt_xy.y);
			A.at<double>(i + 1, 1) = 1 / Z_*(b2*m_fy + b3*pt_xy.y);
			A.at<double>(i + 1, 2) = 1 / Z_*(c2*m_fy + c3*pt_xy.y);
			A.at<double>(i + 1, 3) = -pt_xy.x*sin(m_omega) - (pt_xy.y / m_fy*(pt_xy.x*cos(m_kappa) - pt_xy.y*sin(m_kappa)) - m_fy*sin(m_kappa))*cos(m_omega);
			A.at<double>(i + 1, 4) = -m_fy*cos(m_kappa) - pt_xy.y / m_fy*(pt_xy.x*sin(m_kappa) + pt_xy.y*cos(m_kappa));
			A.at<double>(i + 1, 5) = -pt_xy.x;
			A.at<double>(i + 1, 6) = 0;
			A.at<double>(i + 1, 7) = pt_xy.y / m_fy;
			A.at<double>(i + 1, 8) = 0;
			A.at<double>(i + 1, 9) = 1;
			A.at<double>(i + 1, 10) = pt_xy.y*r2;
			A.at<double>(i + 1, 11) = pt_xy.y*r2*r2;
			A.at<double>(i + 1, 12) = 2 * pt_xy.x*pt_xy.y;
			A.at<double>(i + 1, 13) = r2 + 2 * pt_xy.y*pt_xy.y;

			double x_0 = -m_fx*X_ / Z_ + pt_xy.x*(m_k1*r2 + m_k2*r2*r2) + m_p1*(r2 + 2 * pt_xy.x*pt_xy.x) + 2 * m_p2*pt_xy.x*pt_xy.y + m_x0;
			double y_0 = -m_fy*Y_ / Z_ + pt_xy.y*(m_k1*r2 + m_k2*r2*r2) + m_p2*(r2 + 2 * pt_xy.y*pt_xy.y) + 2 * m_p1*pt_xy.x*pt_xy.y + m_y0;

			L.at<double>(i, 0) = pt_uv.x - x_0;
			L.at<double>(i + 1, 0) = pt_uv.y - y_0;
		}

		cv::Mat_<double> X = (A.t()*A).inv()*(A.t()*L);
		m_X += X.at<double>(0, 0);
		m_Y += X.at<double>(1, 0);
		m_Z += X.at<double>(2, 0);
		m_phi += X.at<double>(3, 0);
		m_omega += X.at<double>(4, 0);
		m_kappa += X.at<double>(5, 0);
		m_fx += X.at<double>(6, 0);
		m_fy += X.at<double>(7, 0);
		m_x0 += X.at<double>(8, 0);
		m_y0 += X.at<double>(9, 0);
		m_k1 += X.at<double>(10, 0);
		m_k2 += X.at<double>(11, 0);
		m_p1 += X.at<double>(12, 0);
		m_p2 += X.at<double>(13, 0);

		++iterate_times;
		
		if (fabs(X.at<double>(3, 0)) < 0.00003
			&&fabs(X.at<double>(4, 0)) < 0.00003
			&&fabs(X.at<double>(5, 0)) < 0.00003)
		{
			m_corrected_elems = { m_fx, m_fy, m_x0, m_y0 ,m_k1, m_k2, m_p1, m_p2 ,m_X, m_Y, m_Z, radianToAngle(m_phi), radianToAngle(m_omega), radianToAngle(m_kappa) };

			double m0 = 0;
			accuracyCompute(m0);

			storage(m0, file);

			std::clog << "共迭代: " << iterate_times - 1 << " 次" << std::endl;
			return m_corrected_elems;
		}

		if (iterate_times > 100)
		{
			std::clog << "收敛失败！" << std::endl;

			system("pause");
			exit(0);
		}
	}
}

void SpaceResection::accuracyCompute(double &m0)const
{
	cv::Mat_<double> A(cv::Size(14, m_imagePoints.size() * 2), 0);
	cv::Mat_<double> L(cv::Size(1, m_imagePoints.size() * 2), 0);

	double a1 = cos(m_phi)*cos(m_kappa) - sin(m_phi)*sin(m_omega)*sin(m_kappa);
	double a2 = -cos(m_phi)*sin(m_kappa) - sin(m_phi)*sin(m_omega)*cos(m_kappa);
	double a3 = -sin(m_phi)*cos(m_omega);
	double b1 = cos(m_omega)*sin(m_kappa);
	double b2 = cos(m_omega)*cos(m_kappa);
	double b3 = -sin(m_omega);
	double c1 = sin(m_phi)*cos(m_kappa) + cos(m_phi)*sin(m_omega)*sin(m_kappa);
	double c2 = -sin(m_phi)*sin(m_kappa) + cos(m_phi)*sin(m_omega)*cos(m_kappa);
	double c3 = cos(m_phi)*cos(m_omega);

	for (size_t i = 0; i < A.rows; i += 2)
	{
		double X_ = a1*(m_controlPoints[i / 2].getX() - m_X) + b1*(m_controlPoints[i / 2].getY() - m_Y) + c1*(m_controlPoints[i / 2].getZ() - m_Z);
		double Y_ = a2*(m_controlPoints[i / 2].getX() - m_X) + b2*(m_controlPoints[i / 2].getY() - m_Y) + c2*(m_controlPoints[i / 2].getZ() - m_Z);
		double Z_ = a3*(m_controlPoints[i / 2].getX() - m_X) + b3*(m_controlPoints[i / 2].getY() - m_Y) + c3*(m_controlPoints[i / 2].getZ() - m_Z);

		Point2D pt_uv = m_imagePoints[i / 2], pt_xy;
		pt_xy.x = pt_uv.x - m_x0;
		pt_xy.y = pt_uv.y - m_y0;

		double r2 = pt_xy.x*pt_xy.x + pt_xy.y*pt_xy.y;
		A.at<double>(i, 0) = 1 / Z_*(a1*m_fx + a3*pt_xy.x);
		A.at<double>(i, 1) = 1 / Z_*(b1*m_fx + b3*pt_xy.x);
		A.at<double>(i, 2) = 1 / Z_*(c1*m_fx + c3*pt_xy.x);
		A.at<double>(i, 3) = pt_xy.y*sin(m_omega) - (pt_xy.x / m_fx*(pt_xy.x*cos(m_kappa) - pt_xy.y*sin(m_kappa)) + m_fx*cos(m_kappa))*cos(m_omega);
		A.at<double>(i, 4) = -m_fx*sin(m_kappa) - pt_xy.x / m_fx*(pt_xy.x*sin(m_kappa) + pt_xy.y*cos(m_kappa));
		A.at<double>(i, 5) = pt_xy.y;
		A.at<double>(i, 6) = pt_xy.x / m_fx;
		A.at<double>(i, 7) = 0;
		A.at<double>(i, 8) = 1;
		A.at<double>(i, 9) = 0;
		A.at<double>(i, 10) = pt_xy.x*r2;
		A.at<double>(i, 11) = pt_xy.x*r2*r2;
		A.at<double>(i, 12) = r2 + 2 * pt_xy.x*pt_xy.x;
		A.at<double>(i, 13) = 2 * pt_xy.x*pt_xy.y;

		A.at<double>(i + 1, 0) = 1 / Z_*(a2*m_fy + a3*pt_xy.y);
		A.at<double>(i + 1, 1) = 1 / Z_*(b2*m_fy + b3*pt_xy.y);
		A.at<double>(i + 1, 2) = 1 / Z_*(c2*m_fy + c3*pt_xy.y);
		A.at<double>(i + 1, 3) = -pt_xy.x*sin(m_omega) - (pt_xy.y / m_fy*(pt_xy.x*cos(m_kappa) - pt_xy.y*sin(m_kappa)) - m_fy*sin(m_kappa))*cos(m_omega);
		A.at<double>(i + 1, 4) = -m_fy*cos(m_kappa) - pt_xy.y / m_fy*(pt_xy.x*sin(m_kappa) + pt_xy.y*cos(m_kappa));
		A.at<double>(i + 1, 5) = -pt_xy.x;
		A.at<double>(i + 1, 6) = 0;
		A.at<double>(i + 1, 7) = pt_xy.y / m_fy;
		A.at<double>(i + 1, 8) = 0;
		A.at<double>(i + 1, 9) = 1;
		A.at<double>(i + 1, 10) = pt_xy.y*r2;
		A.at<double>(i + 1, 11) = pt_xy.y*r2*r2;
		A.at<double>(i + 1, 12) = 2 * pt_xy.x*pt_xy.y;
		A.at<double>(i + 1, 13) = r2 + 2 * pt_xy.y*pt_xy.y;

		double x_0 = -m_fx*X_ / Z_ + pt_xy.x*(m_k1*r2 + m_k2*r2*r2) + m_p1*(r2 + 2 * pt_xy.x*pt_xy.x) + 2 * m_p2*pt_xy.x*pt_xy.y + m_x0;
		double y_0 = -m_fy*Y_ / Z_ + pt_xy.y*(m_k1*r2 + m_k2*r2*r2) + m_p2*(r2 + 2 * pt_xy.y*pt_xy.y) + 2 * m_p1*pt_xy.x*pt_xy.y + m_y0;

		L.at<double>(i, 0) = pt_uv.x - x_0;
		L.at<double>(i + 1, 0) = pt_uv.y - y_0;
	}

	cv::Mat_<double> Q = (A.t()*A).inv();

	cv::Mat_<double> X = Q*(A.t()*L);
	cv::Mat_<double> V = A*X - L;

	double vv = 0;
	for (auto i = 0; i < V.rows; ++i)
	{
		double *data = V.ptr<double>(i);
		for (auto j = 0; j < V.cols; ++j)
			vv += data[j] * data[j];
	}

	m0 = sqrt(vv / (2 * m_imagePoints.size() - 14));
}

void SpaceResection::storage(double &m0, const std::string &file)const
{
	std::fstream fout(file, std::ios::out);
	if (!fout)
		return;

	fout << std::fixed << std::setprecision(6);

	fout << "[Xs,Ys,Zs,phi,omega,kappa]" << std::endl;
	fout <<m_X << "   " << m_Y << "   " << m_Z << "   " << radianToAngle(m_phi) << "   " << radianToAngle(m_omega) << "   " << radianToAngle(m_kappa) << std::endl;
	fout << std::endl;

	fout << "[fx,fy,x0,y0]" << std::endl;
	fout <<m_fx << "   " << m_fy << "   " << m_x0 << "   " << m_y0 << std::endl;
	fout << std::endl;

	fout.unsetf(std::ios::fixed);
	fout << "[k1,k2,p1,p2]" << std::endl;
	fout <<m_k1 << "   " << m_k2 << "   " << m_p1 << "   " << m_p2 << std::endl;
	fout << std::endl;

	fout << std::fixed << std::setprecision(6) << "m0: " << m0 << std::endl;
	fout << std::endl;
}
