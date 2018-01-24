#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd::Identity(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 30;
  
  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.
  
  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  n_x_ = 5;
  n_aug_ = 7;
  lambda_ = 3 - n_aug_;

  //set weights
  weights_ = VectorXd(2 * n_aug_ + 1);
  weights_(0) = lambda_ / (lambda_ + n_aug_);
  for (int i = 1; i < weights_.size(); i++)
  {
	  weights_(i) = 1.0 / (2.0 * (lambda_ + n_aug_));
  }

  R_radar_ = MatrixXd(3, 3);
  R_radar_.fill(0.0);
  R_radar_(0, 0) = std_radr_ * std_radr_;
  R_radar_(1, 1) = std_radphi_ * std_radphi_;
  R_radar_(2, 2) = std_radrd_ * std_radrd_;

  R_lidar_ = MatrixXd(2, 2);
  R_lidar_.fill(0.0);
  R_lidar_(0, 0) = std_laspx_ * std_laspx_;
  R_lidar_(1, 1) = std_laspy_ * std_laspy_;
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
	/**
	TODO:

	Complete this function! Make sure you switch between lidar and radar
	measurements.
	*/
	if (!is_initialized_)
	{
		if (meas_package.sensor_type_ == MeasurementPackage::SensorType::LASER && !use_laser_)
		{
			return;
		}
		if (meas_package.sensor_type_ == MeasurementPackage::SensorType::RADAR && !use_radar_)
		{
			return;
		}

		VectorXd x_aug = VectorXd(n_aug_);
		x_aug.fill(0.0);

		time_us_ = meas_package.timestamp_;
		if (meas_package.sensor_type_ == MeasurementPackage::SensorType::LASER)
		{
			x_(0) = x_aug(0) = meas_package.raw_measurements_(0);
			x_(1) = x_aug(1) = meas_package.raw_measurements_(1);
		}
		else // radar
		{
			double rho = meas_package.raw_measurements_(0);
			double phi = meas_package.raw_measurements_(1);
			x_(0) = x_aug(0) = rho * cos(phi);
			x_(1) = x_aug(1) = rho * sin(phi);
		}

		MatrixXd P_aug_ = MatrixXd(5, 5);
		//create augmented covariance matrix
		P_aug_.fill(0);
		P_aug_.topLeftCorner(n_x_, n_x_) = P_;
		P_aug_(n_x_, n_x_) = std_a_ * std_a_;
		P_aug_(n_x_ + 1, n_x_ + 1) = std_yawdd_ * std_yawdd_;

		//create square root matrix
		MatrixXd A = P_aug_.llt().matrixL();

		//create augmented sigma points
		Xsig_pred_.col(0) = x_aug;
		//set remaining sigma points
		for (int i = 0; i < n_aug_; i++)
		{
			Xsig_pred_.col(i + 1) = x_aug + sqrt(lambda_ + n_aug_) * A.col(i);
			Xsig_pred_.col(i + 1 + n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * A.col(i);
		}

		is_initialized_ = true;
		return;
	}
	double dt = (meas_package.timestamp_ - time_us_) / 1000000.0;
	time_us_ = meas_package.timestamp_;
	Prediction(dt);
	if (meas_package.sensor_type_ == MeasurementPackage::SensorType::LASER)
	{
		UpdateLidar(meas_package);
	}
	else
	{
		UpdateRadar(meas_package);
	}
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
	const double EPS = 0.001;
	VectorXd x_aug = VectorXd(n_aug_);
	P_.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++)
	{
		//predict sigma points
		MatrixXd Xsig_aug = MatrixXd(Xsig_pred_);
		VectorXd x_kn(n_x_);
		VectorXd x_k = Xsig_aug.col(i);
		double p_x = x_k[0];
		double p_y = x_k[1];
		double v = x_k[2];
		double psi = x_k[3];
		double psi_dot = x_k[4];
		double nu_a = x_k[5];
		double nu_psi = x_k[6];

		double nufac = 0.5 * delta_t * delta_t * nu_a;
		if (fabs(psi_dot) > EPS)
		{
			double cosinarg = psi + delta_t * psi_dot;
			x_kn[0] = x_k[0] + (v / psi_dot) * (sin(cosinarg) - sin(psi));
			x_kn[1] = x_k[1] + (v / psi_dot) * (-cos(cosinarg) + cos(psi));
		}
		else //avoid division by zero
		{
			x_kn[0] = x_k[0] + v * cos(psi) * delta_t;
			x_kn[1] = x_k[1] + v * sin(psi) * delta_t;
		}
		x_kn[0] += nufac * cos(psi);
		x_kn[1] += nufac * sin(psi);
		x_kn[2] = x_k[2] + delta_t * nu_a;
		x_kn[3] = x_k[3] + psi_dot * delta_t + 0.5 * delta_t * delta_t * nu_psi;
		x_kn[4] = x_k[4] + delta_t * nu_psi;
		//write predicted sigma points into right column
		Xsig_pred_.col(i) = x_kn;
	}
	// Predict mean and covariance from sigma points
	//predict state mean
	for (int i = 0; i < weights_.size(); i++)
	{
		x_aug += Xsig_pred_.col(i) * weights_(i);
	}
	//predict state covariance matrix
	for (int i = 0; i < weights_.size(); i++)
	{
		VectorXd tmp = Xsig_pred_.col(i) - x_aug;
		while (tmp(3) > M_PI)
			tmp(3) -= 2.*M_PI;
		while (tmp(3) < -M_PI)
			tmp(3) += 2.*M_PI;
		P_ += weights_(i) * tmp * tmp.transpose();
	}
	x_ = x_aug.head(n_x_);
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
	int n_z = 2;
	//create matrix for sigma points in measurement space
	MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

	//mean predicted measurement
	VectorXd z_pred = VectorXd(n_z);
	//measurement covariance matrix S
	MatrixXd S = MatrixXd(n_z, n_z);
	z_pred.fill(0.0);
	S.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++)
	{
		//transform sigma points into measurement space
		double px = Xsig_pred_(0, i);
		double py = Xsig_pred_(1, i);

		Zsig(0, i) = px;
		Zsig(1, i) = py;
	}
	for (int i = 0; i < 2 * n_aug_ + 1; i++)
	{
		//calculate mean predicted measurement
		z_pred += weights_(i) * Zsig.col(i);
	}
	for (int i = 0; i < 2 * n_aug_ + 1; i++)
	{
		//calculate innovation covariance matrix S
		VectorXd z_diff = Zsig.col(i) - z_pred;
		S += weights_(i) * z_diff * z_diff.transpose();
	}
	S += R_radar_;
	//create matrix for cross correlation Tc
	MatrixXd Tc = MatrixXd(n_x_, n_z);

	Tc.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++)
	{
		//calculate cross correlation matrix
		VectorXd x_diff = Xsig_pred_.col(i) - x_;
		VectorXd z_diff = Zsig.col(i) - z_pred;

		Tc += weights_(i) * x_diff * z_diff.transpose();
	}
	//calculate Kalman gain K;
	MatrixXd S_inv = S.inverse();
	MatrixXd K = Tc * S_inv;
	//update state mean and covariance matrix
	VectorXd z_diff = (meas_package.raw_measurements_ - z_pred);
	x_ += K * z_diff;
	P_ -= K * S * K.transpose();

	// TODO: Calculate NIS
	double epsilon = z_diff.transpose() * S_inv * z_diff;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
	int n_z = 3;
	//create matrix for sigma points in measurement space
	MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

	//mean predicted measurement
	VectorXd z_pred = VectorXd(n_z);
	//measurement covariance matrix S
	MatrixXd S = MatrixXd(n_z, n_z);
	z_pred.fill(0.0);
	S.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++)
	{
		//transform sigma points into measurement space
		double px = Xsig_pred_(0, i);
		double py = Xsig_pred_(1, i);
		double v = Xsig_pred_(2, i);
		double psi = Xsig_pred_(3, i);

		double rho = sqrt(px * px + py * py);
		double phi = atan2(py, px);
		double rho_dot = (px * cos(psi) * v + py * sin(psi) * v) / rho;

		Zsig(0, i) = rho;
		Zsig(1, i) = phi;
		Zsig(2, i) = rho_dot;
	}
	for (int i = 0; i < 2 * n_aug_ + 1; i++)
	{
		//calculate mean predicted measurement
		z_pred += weights_(i) * Zsig.col(i);
	}
	for (int i = 0; i < 2 * n_aug_ + 1; i++)
	{
		//calculate innovation covariance matrix S
		VectorXd z_diff = Zsig.col(i) - z_pred;
		while (z_diff(1) > M_PI)
			z_diff(1) -= 2 * M_PI;
		while (z_diff(2) < -M_PI)
			z_diff(1) += 2 * M_PI;
		S += weights_(i) * z_diff * z_diff.transpose();
	}
	S += R_radar_;

	//create matrix for cross correlation Tc
	MatrixXd Tc = MatrixXd(n_x_, n_z);

	Tc.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++)
	{
		//calculate cross correlation matrix
		VectorXd x_diff = Xsig_pred_.col(i) - x_;
		VectorXd z_diff = Zsig.col(i) - z_pred;
		//angle normalization
		while (z_diff(1)> M_PI)
			z_diff(1) -= 2.*M_PI;
		while (z_diff(1)<-M_PI)
			z_diff(1) += 2.*M_PI;

		while (x_diff(3)> M_PI)
			x_diff(3) -= 2.*M_PI;
		while (x_diff(3)<-M_PI)
			x_diff(3) += 2.*M_PI;

		Tc += weights_(i) * x_diff * z_diff.transpose();
	}
	//calculate Kalman gain K;
	MatrixXd S_inv = S.inverse();
	MatrixXd K = Tc * S_inv;
	//update state mean and covariance matrix
	VectorXd z_diff = (meas_package.raw_measurements_ - z_pred);
	//angle normalization
	while (z_diff(1)> M_PI)
		z_diff(1) -= 2.*M_PI;
	while (z_diff(1)<-M_PI)
		z_diff(1) += 2.*M_PI;
	x_ += K * z_diff;
	P_ -= K * S * K.transpose();

	// TODO: Calculate NIS and save NIS
	double epsilon = z_diff.transpose() * S_inv * z_diff;
}
