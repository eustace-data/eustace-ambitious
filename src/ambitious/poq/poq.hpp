#ifndef POQ_HPP
#define POQ_HPP

#include "ambitious/versions/versions.hpp"
VERSIONS_ADD(poq_hpp, "$Revision: 1235 $")

#include <map>
#include <Eigen/Dense>
// #include <iostream>
#include <boost/random.hpp>
#include <boost/math/distributions/normal.hpp>


template<typename Scalar>
Scalar logit(Scalar p) {
  return std::log(p) - std::log(1-p);
}
template<typename Scalar>
Scalar ilogit(Scalar x) {
  return 1.0/(1.0 + std::exp(-x));
}

template <typename DerivedA, typename DerivedB>
void logit(const Eigen::DenseBase<DerivedA>& p,
	   Eigen::DenseBase<DerivedB> const & x_) {
  Eigen::DenseBase<DerivedB>& x = const_cast< Eigen::DenseBase<DerivedB>& >(x_);
  x = p.derived().array().log() - (1.0 - p.derived().array()).log();
}
template <typename DerivedA, typename DerivedB>
void ilogit(const Eigen::DenseBase<DerivedA>& x,
	   Eigen::DenseBase<DerivedB> const & p_) {
  Eigen::DenseBase<DerivedB>& p = const_cast< Eigen::DenseBase<DerivedB>& >(p_);
  p = 1.0/(1.0 + (-x.derived().array()).exp());
}








template<typename Scalar_>
class poq_distribution {
public:
  typedef Scalar_ Scalar;
  typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Vector;
  typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> Matrix;
  typedef Vector param_type;
protected:
  /* limits for direct vs Taylor series noise vs bias tradeoff */
#ifdef OLDCXX
#define STATIC_CONSTEXPR_COMPATIBILITY
#else
#define STATIC_CONSTEXPR_COMPATIBILITY static constexpr
#endif
  /* Approx 2.961809e-4 for double. Limit for theta. */
  STATIC_CONSTEXPR_COMPATIBILITY Scalar taylor_limit_ =
    std::exp(std::log(std::numeric_limits<Scalar>::epsilon() * 8.0) / 4.0) / std::log(2.0);
  /* Approx 2.701860e-4 for double. Limit for theta*log(2*p). */
  STATIC_CONSTEXPR_COMPATIBILITY Scalar taylor_limit_forward_ = 
    std::exp(std::log(std::numeric_limits<Scalar>::epsilon() * 24.0) / 4.0);
  /* Approx 9.612435e-6 for double. Limit for 2*theta*q. */
  STATIC_CONSTEXPR_COMPATIBILITY Scalar taylor_limit_inverse_ =
    std::exp(std::log(std::numeric_limits<Scalar>::epsilon() * 4.0) / 3.0);

  typename boost::random::uniform_real_distribution<Scalar> unif_;
  param_type theta_;
  param_type theta_real_;
  Eigen::Matrix<Scalar, 2, 1> limits_minus_;
  Eigen::Matrix<Scalar, 2, 1> limits_plus_;
  Eigen::Matrix<bool, 2, 1> limits_finite_;
  Eigen::Matrix<Scalar, 2, 1> limits_;

  Scalar tailMinus(Scalar p, Scalar theta) const {
    Scalar q;
    // Compute tail q = f^-(p)
    Scalar u = std::log(2.0 * p);
    if (std::fabs(theta * u) < taylor_limit_forward_) {
      q = (u / 2.0) * (1.0 - (theta/2.0) * u * (1.0 - (theta/3.0) * u));
    } else {
      q = (1.0 - std::exp(-theta * u)) / (2.0 * theta);
    }
    return q;
  }

  Scalar tailPlus(Scalar p, Scalar theta) const {
    return -tailMinus(1.0-p, theta);
  }

  Scalar tailMinusInverse(Scalar q, Scalar theta) const {
    Scalar p;
    // Compute tail inverse p, such that q = f^-(p)
    Scalar v = 2.0 * theta * q;
    if (std::fabs(v) < taylor_limit_inverse_) {
      p = std::exp(q * (2.0 + v * (1.0 + v * (2.0 / 3.0)))) / 2.0;
    } else if (q >= tailMinus(1.0, theta)) { // Upper limit reached
      p = 1.0;
    } else if ((theta > 0) | ((theta < 0) & (v < 1.0))) { // No lower limit, or not reached
      p = std::exp(std::log(1.0 - v) / (-theta)) / 2.0;
    } else {
      p = 0.0;
    }
    return p;
  }

  Scalar tailPlusInverse(Scalar q, Scalar theta) const {
    return 1.0-tailMinusInverse(-q, theta);
  }

  template <typename Derived, typename OtherDerived>
  void tailMinus(const Eigen::DenseBase<Derived>& p,
		 Scalar theta,
		 Eigen::DenseBase<OtherDerived> const & f_) const {
    Eigen::DenseBase<OtherDerived>& f = const_cast< Eigen::DenseBase<OtherDerived>& >(f_);
    f.derived().resize(p.rows(), p.cols());
    for (int i=0; i < p.size(); ++i) {
      f.derived()[i] = tailMinus(p.derived()[i], theta);
    }
  }

  template <typename Derived, typename OtherDerived>
  void tailMinusInverse(const Eigen::DenseBase<Derived>& q,
			Scalar theta,
			Eigen::DenseBase<OtherDerived> const & p_) const {
    Eigen::DenseBase<OtherDerived>& p = const_cast< Eigen::DenseBase<OtherDerived>& >(p_);
    p.derived().resize(q.rows(), q.cols());
    for (int i=0; i < q.size(); ++i) {
      p.derived()[i] = tailMinusInverse(q.derived()[i], theta);
    }
  }
  template <typename Derived, typename OtherDerived>
  void tailPlus(const Eigen::DenseBase<Derived>& p,
		Scalar theta,
		Eigen::DenseBase<OtherDerived> const & f_) const {
    Eigen::DenseBase<OtherDerived>& f = const_cast< Eigen::DenseBase<OtherDerived>& >(f_);
    f.derived().resize(p.rows(), p.cols());
    for (int i=0; i < p.size(); ++i) {
      f.derived()[i] = tailPlus(p.derived()[i], theta);
    }
  }
  
  template <typename Derived, typename OtherDerived>
  void tailPlusInverse(const Eigen::DenseBase<Derived>& q,
		       Scalar theta,
		       Eigen::DenseBase<OtherDerived> const & p_) const {
    Eigen::DenseBase<OtherDerived>& p = const_cast< Eigen::DenseBase<OtherDerived>& >(p_);
    p.derived().resize(q.rows(), q.cols());
    for (int i=0; i < q.size(); ++i) {
      p.derived()[i] = tailPlusInverse(q.derived()[i], theta);
    }
  }

  template <typename DerivedA, typename DerivedB>
  void tailMinusDeriv(const Eigen::DenseBase<DerivedA>& p,
		      Scalar theta,
		      int Dp, int Dt, bool Drho,
		      Eigen::DenseBase<DerivedB> const & f_) const {
    Eigen::DenseBase<DerivedB>& f = const_cast< Eigen::DenseBase<DerivedB>& >(f_);
    //    f.derived().resize(p.rows(), p.cols());

    if (Dt == 0) {
      if (Dp == 0) {
	// Compute tail f^-(p)
	Vector tmp;
	tailMinus(p, theta, tmp);
	f = tmp;
      } else if (Dp == 1) {
	f = (2.0*p.derived().array()).pow(-theta-1);
	if (Drho) {
	  f = f.derived().cwiseProduct(p.derived().cwiseProduct((1.0-p.derived().array()).matrix()));
	}
      } else { // Dp >= 2, Dt == 0
	if (!Drho) {
	  Scalar thPower = theta + 1.0;
	  for (int k=1; k <= Dp-2; ++k) {
	    thPower *= theta + k + 1;
	  }
	  f = (std::log(2.0)*(-theta-1) + (-theta-Dp)*p.derived().array().log()).exp() * thPower;
	  if (Dp % 2 == 0) {
	    f = -f.derived();
	  }
	} else { // Drho == true
	  if (Dp > 2) {
	    std::cout << "Dp = " << Dp << ", Dt = " << Dt << std::endl;
	  }
	  assert(Dp <= 2);
	  // Dp == 2
	  // log(2p)
	  Vector logP = (2.0*p.derived().array()).log().matrix();
	  // p(1-p)
	  Vector PP = p.derived().cwiseProduct((1.0-p.derived().array()).matrix());
	  // (2p)^(-theta-1)
	  Matrix tmp = ((-theta-1)*logP.array()).exp().matrix();
	  // f = -2(theta+1) (2p)^(-theta-Dp)
	  f = (-2*(theta+1)) * (logP*(-theta-Dp)).array().exp();
	  // f = (f p (1-p) + (2p)^(-theta-1) (1-2p)) p (1-p)
	  f = (f.derived().cwiseProduct(PP) +
	       tmp.cwiseProduct((1-2.0*p.derived().array()).matrix())).cwiseProduct(PP);
	}
      }
    } else if (Dp == 0) { // Dt > 0
      Vector logP = (2.0*p.derived().array()).log();
      if (Dt == 1) {
	if (std::fabs(theta) < taylor_limit_) {
	  // thLog = theta * pmax(-700, log(2 * p)); // -700 approx log(1e-304)
	  f = (logP.cwiseProduct(logP).array() / (-4.0)) *
	    (1.0 - (2*theta/3.0) * logP.array()).array();
	} else {
	  f = (-theta * logP).array().exp() * (1.0 + (theta * logP).array() - 1.0) / (2*theta*theta);
	}
      } else if (Dt == 2) {
	if (std::fabs(theta) < taylor_limit_) {
	  f = (3.0*logP.array()).exp() / 6.0;
	} else {
	  Matrix tmp = (1.0 + (theta * logP).array()).matrix();
	  f = ((-theta * logP).array().exp() * (1 + tmp.cwiseProduct(tmp).array()) - 2.0) / (-2*theta*theta*theta);
	}
      } else {
	std::cout << "Dp = " << Dp << ", Dt = " << Dt << " > 2" << std::endl;
	assert(false);
      }
    } else { // Dp > 0, Dt > 0
      Vector logP = (2.0*p.derived().array()).log();
      Matrix w1 = logP;
      for (int k=1; k < Dp; ++k) {
	w1.array() += 1/(theta+k);
      }
      Matrix ff;
      tailMinusDeriv(p, theta, Dp, 0, false, ff);
      if (Dt == 1) {
	f = w1.cwiseProduct(ff);
      } else if (Dt == 2) {
	Scalar w2 = 0.0;
	for (int k=1; k < Dp; ++k) {
	  w2 += 1/(theta+k)/(theta+k);
	}
	f = (w1.cwiseProduct(w1).array() - w2) * ff.array();
      } else {
	std::cout << "Dp = " << Dp << ", Dt = " << Dt << " > 2" << std::endl;
	assert(false);
      }
      if (Drho) {
	// p(1-p)
	Vector PP = p.derived().cwiseProduct((1.0-p.derived().array()).matrix());
	if (Dt == 1) {
	  f = f.derived().cwiseProduct(PP);
	} else if (Dt == 2) {
	  tailMinusDeriv(p, theta, 1, 1, false, ff);
	  // (2p)^(-theta-1)
	  // f = (f p (1-p) + ff (1-2p)) p (1-p)
	  f = (f.derived().cwiseProduct(PP) +
	       ff.cwiseProduct((1.0-2.0*p.derived().array()).matrix())).cwiseProduct(PP);
	} else {
	  std::cout << "Dp = " << Dp << ", Dt = " << Dt << " > 2" << std::endl;
	  assert(false);
	}
      }
    }
  }

  template <typename DerivedA, typename DerivedB>
  void tailPlusDeriv(const Eigen::DenseBase<DerivedA>& p,
		      Scalar theta,
		      int Dp, int Dt, bool Drho,
		      Eigen::DenseBase<DerivedB> const & f_) const {
    Eigen::DenseBase<DerivedB>& f = const_cast< Eigen::DenseBase<DerivedB>& >(f_);
    //    f.derived().resize(p.rows(), p.cols());

    // Compute tail Df^+(p)/D...=(-)^(Dp+1) Df^-(1-p)/D...
    tailMinusDeriv((1.0-p.derived().array()).matrix(), theta, Dp, Dt, Drho, f);
    if (Dp % 2 == 0) {
      f = -f.derived();
    }
  }


  void calc_limits() {
    for (int i=0; i < 2; ++i) {
      limits_finite_[i] = limitIsFinite(theta_real_(3+i));
      limits_minus_[i] = limitMinus(i==1, theta_real_(3));
      limits_plus_[i] = limitPlus(i==1, theta_real_(4));
      limits_[i] = theta_real_(0) + theta_real_(1) * ((1.0-theta_real_(2)) * limits_minus_[i] +
						    theta_real_(2) * limits_plus_[i]);
    }
    if (limits_finite_[0]) {
      limits_[0] = theta_real_(0) + theta_real_(1) * ((1.0-theta_real_(2)) * limits_minus_[0] +
						    theta_real_(2) * limits_plus_[0]);
    } else {
      limits_[0] = -std::numeric_limits<Scalar>::infinity();
    }
    if (limits_finite_[1]) {
      limits_[1] = theta_real_(0) + theta_real_(1) * ((1.0-theta_real_(2)) * limits_minus_[1] +
						    theta_real_(2) * limits_plus_[1]);
    } else {
      limits_[1] = std::numeric_limits<Scalar>::infinity();
    }
  }

  void calc_real_theta() {
    theta_real_ = theta_;
    theta_real_(1) = exp(theta_(1));
    theta_real_(2) = ilogit(theta_(2));
    calc_limits();
  };

public:
  poq_distribution() : unif_(0.0, 1.0), theta_(param_type::Zero(5)) {
    calc_real_theta();
  };
  poq_distribution(const param_type& theta) : unif_(0.0, 1.0), theta_(theta) {
    calc_real_theta();
  };
  void param(const param_type& theta) {
    theta_ = theta;
    calc_real_theta();
  };
  const param_type& param() const {
    return theta_;
  };
  const param_type& param_real() const {
    return theta_real_;
  };



  bool limitIsFinite(Scalar theta) {
    return (theta < 0.0);
  }
  Scalar limitMinus(bool upper, Scalar theta) {
    if (upper) {
      return tailMinus(1.0, theta);
    } else if (limitIsFinite(theta)) {
      return 1.0/(2.0*theta);
    } else {
      return -std::numeric_limits<Scalar>::infinity();
    }
  }
  Scalar limitPlus(bool upper, Scalar theta) {
    if (!upper) {
      return tailPlus(0.0, theta);
    } else if (limitIsFinite(theta)) {
      return -1.0/(2.0*theta);
    } else {
      return std::numeric_limits<Scalar>::infinity();
    }
  }
  bool limitIsFinite(int upper) const {
    return limits_finite_[upper];
  }
  Scalar limitMinus(int upper) const {
    return limits_minus_[upper];
  }
  Scalar limitPlus(int upper) const {
    return limits_plus_[upper];
  }
  Scalar limit(int upper) const {
    return limits_[upper];
  }
  Scalar taylorLimit() const {
    return taylor_limit_;
  }
  Scalar taylorLimitForward() const {
    return taylor_limit_forward_;
  }
  Scalar taylorLimitInverse() const {
    return taylor_limit_inverse_;
  }

/** 
 * 
 * 
 * @param p 
 * @param q 
 */
  template <typename Derived, typename OtherDerived>
  void quantile(const Eigen::DenseBase<Derived>& p,
		Eigen::DenseBase<OtherDerived> const & q_) const {
    Eigen::DenseBase<OtherDerived>& q = const_cast< Eigen::DenseBase<OtherDerived>& >(q_);
    //    q.derived().resize(p.rows());
    
    Scalar gammaPlus = theta_real_(2);
    Scalar gammaMinus = 1-theta_real_(2);
    Derived fminus, fplus;
    tailMinus(p, theta_real_(3), fminus);
    tailPlus(p, theta_real_(4), fplus);
    q.derived().array() = theta_real_(0) + ((theta_real_(1) * gammaMinus) * fminus.derived() +
			 (theta_real_(1) * gammaPlus) * fplus.derived()).array();
  }
  Scalar quantile(Scalar p) const {
    Scalar fminus = tailMinus(p, theta_real_(3));
    Scalar fplus = tailPlus(p, theta_real_(4));
    return theta_real_(0) + theta_real_(1) * ((1.0 - theta_real_(2)) * fminus + theta_real_(2) * fplus);
  }


  /** 
   * 
   * 
   * @param p 
   * @param D 
   * @param Dp 
   * @param Dt 
   */
  template <typename DerivedA, typename DerivedB, typename DerivedC>
  void quantileDerivative(const Eigen::DenseBase<DerivedA>& p,
			  Eigen::DenseBase<DerivedB> const & D_,
			  int Dp, int Dt, bool Drho,
			  const Eigen::DenseBase<DerivedC>& reduce) const  {
    Eigen::DenseBase<DerivedB>& D = const_cast< Eigen::DenseBase<DerivedB>& >(D_);
    D.derived().resize(p.rows(), p.cols());

    Scalar gammaPlus = theta_real_(2);
    Scalar gammaMinus = 1.0-theta_real_(2);
    if (Dt == 0) {
      Vector Dminus0, Dplus0;
      tailMinusDeriv(p, theta_real_(3), Dp, 0, Drho, Dminus0);
      tailPlusDeriv(p, theta_real_(4), Dp, 0, Drho, Dplus0);
      if (reduce.size() <= 1) {
	D = (theta_real_(1) * gammaMinus) * Dminus0;
	D += (theta_real_(1) * gammaPlus) * Dplus0;
	D.derived().array() += theta_real_(0) * (Dp == 0);
      } else { // reduce is non-empty
	D.derived().resize(1, 1);
	D.derived()(0,0) = (theta_real_(1) * gammaMinus) * Dminus0.dot(reduce.derived());
	D.derived()(0,0) += (theta_real_(1) * gammaPlus) * Dplus0.dot(reduce.derived());
	D.derived()(0,0) += theta_real_(0) * (Dp == 0);
      }
    } else if (Dt == 1) {
      Vector Dminus0, Dplus0, Dminus1, Dplus1;
      tailMinusDeriv(p, theta_real_(3), Dp, 0, Drho, Dminus0);
      tailMinusDeriv(p, theta_real_(3), Dp, 1, Drho, Dminus1);
      tailPlusDeriv(p, theta_real_(4), Dp, 0, Drho, Dplus0);
      tailPlusDeriv(p, theta_real_(4), Dp, 1, Drho, Dplus1);
      if (reduce.size() <= 1) {
	D.derived().resize(p.rows(), 5);
	D.col(0).array() = 1.0*(Dp == 0);
	D.col(1) = (theta_real_(1) * gammaMinus) * Dminus0;
	D.col(1) += (theta_real_(1) * gammaPlus) * Dminus0;
	D.col(2) = (theta_real_(1) * gammaMinus * gammaPlus) *
	  (Dplus0 - Dminus0);
	D.col(3) = (theta_real_(1) * gammaMinus) * Dminus1;
	D.col(4) = (theta_real_(1) * gammaPlus) * Dminus1;
      } else { // reduce != NULL
	D.derived().resize(1, 5);
	D.derived()(0,0) = 1.0*(Dp == 0);
	D.derived()(0,1) = (theta_real_(1) * gammaMinus) * Dminus0.dot(reduce.derived());
	D.derived()(0,1) += (theta_real_(1) * gammaPlus) * Dminus0.dot(reduce.derived());
	D.derived()(0,2) = (theta_real_(1) * gammaMinus * gammaPlus) *
	  (Dplus0.dot(reduce.derived()) - Dminus0.dot(reduce.derived()));
	D.derived()(0,3) = (theta_real_(1) * gammaMinus) * Dminus1.dot(reduce.derived());
	D.derived()(0,4) = (theta_real_(1) * gammaPlus) * Dminus1.dot(reduce.derived());
      }
    } else if (Dt == 2) {
      assert(reduce.size() == p.size());
      Scalar Dminus[3], Dplus[3];
      Vector tmp;
      for (int k=0; k < 3; ++k) {
	tailMinusDeriv(p, theta_real_(3), Dp, k, Drho, tmp);
	Dminus[k] = tmp.dot(reduce.derived());
	tailPlusDeriv(p, theta_real_(4), Dp, k, Drho, tmp);
	Dplus[k] = tmp.dot(reduce.derived());
      }
      Scalar bgMinus = theta_real_[1] * gammaMinus;
      Scalar bgPlus = theta_real_[1] * gammaPlus;
      D = Matrix::Zero(5,5);
      D.derived()(1,1) =
	bgMinus * Dminus[0] +
	bgPlus * Dminus[0];
      D.derived()(1,2) = D(2,1) =
	bgMinus * (-gammaPlus * Dminus[0]) +
	bgPlus * (gammaMinus * Dplus[0]);
      D.derived()(1,3) = D(3,1) =
	bgMinus * (Dminus[1]);
      D.derived()(1,4) = D(4,1) =
	bgPlus * (Dplus[1]);
      D.derived()(2,2) = D(2,2) =
	(gammaMinus - gammaPlus) * (
	bgMinus * (-gammaPlus * Dminus[0]) +
	bgPlus * (gammaMinus * Dplus[0]) );
      D(2,3) = D(3,2) =
	bgMinus * (-gammaPlus * Dminus[1]);
      D(2,4) = D(4,2) =
	bgPlus * (gammaMinus * Dplus[1]);
      D(3,3) =
	bgMinus * (Dminus[2]);
      D(4,4) =
	bgPlus * (Dplus[2]);
    } else {
      assert(false);
    }
  }


  template <typename DerivedB>
  void quantileDerivative(Scalar p,
			  Eigen::DenseBase<DerivedB> const & D_,
			  int Dp, int Dt, bool Drho) const {
    Eigen::Matrix<Scalar, 1, 1> p_;
    p_[0] = p;
    quantileDerivative(p_, D_, Dp, Dt, Drho, Eigen::Matrix<Scalar, 1, 1>::Ones());
  }
  
  Scalar quantileDerivativeDp(Scalar p, int Dp, bool Drho) const {
    typename Eigen::Matrix<Scalar, 1, 1> p_;
    typename Eigen::Matrix<Scalar, 1, 1> D_;
    p_[0] = p;
    quantileDerivative(p_, D_, Dp, 0, Drho, Eigen::Matrix<Scalar, 1, 1>::Ones());
    return D_[0];
  }


  Scalar quantileInverse(Scalar q) const {
    /* Find p such that q=quantile(p)
     * Assume that both tail parameters are > -1
     * For q > theta_real_(0),
     * - the solution p1 to q=gammaPlus*tailPlus(p1,tailparameter)
     *   fulfils p < p1, and
     * - the solution p2 to q=gammaMinus*(p2-0.5) + gammaPlus*tailPlus(p2,tailparameter)
     *   fulfils p2 < p.
     * - The 2nd order derivative of tailPlus is positive 
     */

    const Scalar tolerance = std::numeric_limits<Scalar>::epsilon() * 1.0e2;

    // Initial estimate for q > median (and reverse for q < median):
    // 1) Find beta such that beta * tailPlus has the same upper bound as the full quantile function
    // 2) Invert q = beta * tailPlus(p0) to get p0

    //std::cout << "theta = " << theta_.transpose() << std::endl;
    //std::cout << "limits = " << limits_.transpose() << std::endl;
    //std::cout << "cdf(" << q << ") = ";

    Scalar qq = (q - theta_real_(0))/theta_real_(1);
    Scalar p;
    Scalar beta;
    if (qq < 0.0) { // Lower tail
      if (limitIsFinite(0)) {
	if (q <= limit(0)) {
	  //std::cout << 0.0 << std::endl;
	  return 0.0;
	}
	beta = (1.0-theta_real_(2)) + theta_real_(2) * limitPlus(0) / limitMinus(0);
      } else {
	beta = 1.0-theta_real_(2);
      }
      p = tailMinusInverse(qq/beta, theta_real_(3));
    } else if (qq > 0.0) { // Upper tail
      if (limitIsFinite(1)) {
	if (q >= limit(1)) {
	  //std::cout << 1.0 << std::endl;
	  return 1.0;
	}
	beta = (1-theta_real_(2)) * limitMinus(1) / limitPlus(1) + theta_real_(2);
      } else {
	beta = theta_real_(2);
      }
      p = tailPlusInverse(qq/beta, theta_real_(4));
    } else { // At the median!
      //std::cout << 0.5 << std::endl;
      return 0.5;
    }
    // std::cout << "limitIsFinite: " << limitIsFinite(0) << ", " << limitIsFinite(1) << std::endl;
    // std::cout << "Minus limits: " << limitMinus(0) << ", " << limitMinus(1) << std::endl;
    // std::cout << "Plus  limits: " << limitPlus(0) << ", " << limitPlus(1) << std::endl;
    // std::cout << "limit: " << limit(0) << ", " << limit(1) << std::endl;
    // std::cout << "theta_real: " << theta_real_.transpose() << std::endl;
    // std::cout << "q, qq, beta: " << q << ", " << qq << ", " << beta << std::endl;

    // Iterative linear zero seeker with partial quadratic adjustment.
    Scalar residual;
    Scalar D1;
    Scalar D2;
    Scalar delta;
    Scalar D1_1;

    //std::cout << "Finding p for q = " << q << std::endl;
    //std::cout << "p[ ] = " << p;
    const int max_loop = 10;
    for (int loop=0 ; loop < max_loop; ++loop) {
      residual = quantile(p) - q;
      if (std::fabs(residual) < tolerance) {
	break;
      }
      //std::cout << ",\t q = " << quantile(p)
      //	<< ",\t Residual = " << residual << std::endl;
      D1 = quantileDerivativeDp(p, 1, false);
      D2 = quantileDerivativeDp(p, 2, false);
      delta = residual / D1;
      /*
       * residual < 0: Moving up, and quadratic adjustment ok if D2>0 or (D2<0 and local max > q)
       * residual > 0: Moving down, and quadratic adjustment ok if D2<0 or (D2>0 and local min < q)
       * q(p) = q(p0) + D1*(p-p0) + D2*(p-p0)^2/2
       * Local optimum when  D1 + D2*(p-p0) = 0, i.e. p-p0 = -D1/D2, and then
       * q(p0-D1/D2) = q(p0) - D1^2/D2 + D2*D1^2/D2^2/2 = q(p0) - 1/2 * D1^2/D2
       * D2 < 0:  q(p0-D1/D2) > q when
       *   q(p0) - 1/2 * D1^2/D2 > q,
       *   residual - 1/2 * D1^2/D2 > 0,
       *   residual*D2 - 1/2 * D1^2 < 0
       * D2 > 0:  q(p0-D1/D2) < q when
       *   q(p0) - 1/2 * D1^2/D2 < q,
       *   residual - 1/2 * D1^2/D2 < 0,
       *   residual*D2 - 1/2 * D1^2 < 0
       */
      //std::cout << "d = " << delta
      //	<< ",\tq(p-d) = " << quantile(p-delta)
      //	<< ",\tr = " << quantile(p-delta) - q << std::endl;
      if (((residual < 0.0) & ((D2 > 0.0) | (residual * D2 < 0.5 * D1 * D1))) |
	  ((residual > 0.0) & ((D2 < 0.0) | (residual * D2 < 0.5 * D1 * D1)))) {
	D1_1 = D1 - D2 * delta;
	delta += D2 * delta * delta / 2.0 / D1_1;
	//std::cout << "d = " << delta
	//	  << ",\tq(p-d) = " << quantile(p-delta)
	//	  << ",\tr = " << quantile(p-delta) - q << std::endl;
      } else {
	//std::cout << "r = " << residual
	//	  << ",\tD2 = " << D2(0)
	//	  << ",\trD2 = " << residual * D2(0)
	//	  << ",\tD1^2/2 = " << 0.5 * D1(0) * D1(0) << std::endl;
      }
      while ((p <= delta) | (p - 1.0 >= delta)) {
	// std::cout << "delta/2" << std::endl;
	delta /= 2.0;
      }
      p -= delta;
      //std::cout << "p[" << loop << "] = " << p;
    }
    //std::cout << ",\t q = " << quantile(p)
    //	      << ",\t Residual = " << quantile(p) - q << std::endl;

    return p;
  }
  

  

  Scalar cdf(Scalar x) const {
    return quantileInverse(x);
  }

  template <typename DerivedA, typename DerivedB>
  void cdf(const Eigen::DenseBase<DerivedA>& x,
	   Eigen::DenseBase<DerivedB> const & p_) const {
    typename Eigen::DenseBase<DerivedB>& p =
      const_cast< typename Eigen::DenseBase<DerivedB>& >(p_);
    p.derived().resize(x.rows(), x.cols());
    for (int i=0; i < x.size(); ++i) {
      p.derived()[i] = quantileInverse(x.derived()[i]);
    }
  }

  template <bool log=false, bool quantile_input=false>
  Scalar pdf(Scalar x) const {
    if (quantile_input &
	((limitIsFinite(0) & (x <= limit(0))) |
	 (limitIsFinite(1) & (x >= limit(1))))) {
      if (log) {
	return -std::numeric_limits<Scalar>::infinity();
      } else {
	return 0.0;
      }
      return 0.0;
    } else {	    
      Vector tmp(1);
      if (quantile_input) {
	quantileDerivative(cdf(x), tmp, 1, 0, false);
      } else {
	quantileDerivative(x, tmp, 1, 0, false);
      }
      if (log) {
	return -std::log(tmp[0]);
      } else {
	return 1.0 / tmp[0];
      }
    }
  }

  template <bool log=false, bool quantile_input=true,
	    typename DerivedA, typename DerivedB>
  void pdf(const Eigen::DenseBase<DerivedA>& x,
	   Eigen::DenseBase<DerivedB> const & y_) const {
    typename Eigen::DenseBase<DerivedB>& y =
      const_cast< typename Eigen::DenseBase<DerivedB>& >(y_);
    Vector tmp(x.size());
    if (quantile_input) {
      cdf(x, y);
      quantileDerivative(y, tmp, 1, 0, false, Eigen::Matrix<Scalar, 1, 1>::Ones());
    } else {
      y.derived().resize(x.rows(),x.cols());
      quantileDerivative(x, tmp, 1, 0, false, Eigen::Matrix<Scalar, 1, 1>::Ones());
    }
    for (int i=0; i < x.size(); ++i) {
      if (quantile_input &
	  ((limitIsFinite(0) & (x.derived()[i] <= limit(0))) |
	   (limitIsFinite(1) & (x.derived()[i] >= limit(1))))) {
	if (log) {
	  y.derived()[i] = -std::numeric_limits<Scalar>::infinity();
	} else {
	  y.derived()[i] = 0.0;
	}
      } else {
	if (log) {
	  y.derived()[i] = -std::log(tmp[i]);
	} else {
	  y.derived()[i] = 1 / tmp[i];
	}
      }
    }
  }




  template <typename DerivedA, typename DerivedB>
  Scalar scoreLog(const Eigen::DenseBase<DerivedA>& x,
		  Eigen::DenseBase<DerivedB> const & score_) const {
    typename Eigen::DenseBase<DerivedB>& score =
      const_cast< typename Eigen::DenseBase<DerivedB>& >(score_);
    pdf<true>(x, score);
    return score.mean();
  }
  template <typename DerivedA, typename DerivedB>
  Scalar scoreLog(const Eigen::DenseBase<DerivedA>& x) const {
    Vector score;
    return scoreLog(x, score);
  }

  template <typename DerivedA, typename DerivedB, typename DerivedC>
  Scalar scoreLogDerivatives(const Eigen::DenseBase<DerivedA>& x,
			     Eigen::DenseBase<DerivedB> const & scoreDt1_,
			     Eigen::DenseBase<DerivedC> const & scoreDt2_) const {
    typename Eigen::DenseBase<DerivedB>& scoreDt1 =
      const_cast< typename Eigen::DenseBase<DerivedB>& >(scoreDt1_);
    typename Eigen::DenseBase<DerivedC>& scoreDt2 =
      const_cast< typename Eigen::DenseBase<DerivedC>& >(scoreDt2_);

    Vector p;
    Vector tmp;

    cdf(x, p);
    pdf<true,false>(p, tmp);

    /* From R code:
     *
     * DS <- poq.score.log(theta=theta, p=p, x.range=x.range,
     *                     derivatives=TRUE, rho=TRUE)
     * Dr <- rho.derivatives(theta, p, collectDtt=DS$Dp)
     * Dt <- DS$Dt + colSums(DS$Dp * Dr$Dt)
     * tmp <- t(DS$Dpt) %*% Dr$Dt
     * Dtt <- DS$Dtt + Dr$Dtt + tmp + t(tmp) + t(Dr$Dt) %*% (DS$Dpp * Dr$Dt)
     * list(value=DS$value, Dt=Dt, Dtt=Dtt)
     */

    /*
     *  Dt1 = dS/d\theta^\top
     *      = \partial{S}/\partial\theta^\top
     *        + \partial{S}/\partial{p}  dp/d\theta^\top
     *      = dSt + dSp dPt
     *  Dt2 = d^2{S}/{d\theta d\theta^\top}
     *      = \partial^2{S}/{\partial\theta\partial\theta^\top}
     *        + \partial^2{S}/{\partial{p}^\top\partial\theta}  dp/d\theta^\top
     *        + dp^\top/d\theta  \partial^2{S}/{\partial{p}\partial\theta^\top}
     *        + dp^\top/d\theta \partial^2{S}/{\partial{p}\partial{p}^\top} dp/d\theta^\top
     *        + \sum_{k=1}^n (\partial{S}/\partial p)_k
     *                       \partial^2{p_k}/{\partial\theta\partial\theta^\top}
     *      = dStt + dSpt dPt + dPt DSpt + dPt dSpp dPt + reduce(dPtt, dSp)
     * Required:
     *   dSt, dStt, dSp, dSpt, dSpp, dPt, reduce(dPtt, dSp)
     */

    
    //    score = tmp;
    return tmp.mean();
  }



  /*
   * Methods for compatibility with boost::random distribution classes
   */

  Scalar mean() const {
    if (theta_real_[3] >= 1.0) {
      if (theta_real_[4] < 1.0) {
	return -std::numeric_limits<Scalar>::infinity();
      } else {
	return std::numeric_limits<Scalar>::quiet_NaN();
      }
    } else if (theta_real_[4] >= 1.0) {
      return std::numeric_limits<Scalar>::infinity();
    }
    return theta_real_[0] + theta_real_[1] *
      ((1.0-theta_real_[2]) / (1.0 - theta_real_[3]) *
       (tailMinus(1) - 0.5)
       +
       theta_real_[2] / (1.0 - theta_real_[4]) *
       (tailPlus(0) + 0.5)
       );
  }
  Scalar sigma() const {
    if ((theta_real_[3] >= 0.5) | (theta_real_[4] >= 0.5)) {
      return std::numeric_limits<Scalar>::infinity();
    }
    return theta_real_[1] *
      ((1.0-theta_real_[2]) / (1.0 - theta_real_[3]) *
       (tailMinus(1) - 0.5)
       +
       theta_real_[2] / (1.0 - theta_real_[4]) *
       (tailPlus(0) + 0.5)
       );
  }
  Scalar min() const {
    return limits_[0];
  }
  Scalar max() const {
    return limits_[1];
  }
  void reset() {
    unif_.reset();
  }
  
  template <class Engine>
  Scalar operator()(Engine& eng) {
    return quantile(unif_(eng));
  }

  template <class Engine>
  Scalar operator()(Engine& eng, const param_type& param) {
    return poq_distribution<Scalar>(param)(eng);
  }


};








/* Vectorise boost::math distribution functions */
template <class Dist, class DerivedA, class DerivedB>
void quantile(const Dist& dist,
	      const Eigen::DenseBase<DerivedA>& x,
	      Eigen::DenseBase<DerivedB> const & y_) {
  typename Eigen::DenseBase<DerivedB>& y =
    const_cast< typename Eigen::DenseBase<DerivedB>& >(y_);
  y.derived().resize(x.rows(), x.cols());
  for (int i=0; i < x.size(); ++i) {
    //        std::cout << "quantile(dist, " << x(i) << ")";
    y.derived()[i] = quantile(dist, x.derived()[i]);
    //        std::cout << " = " << y(i) << std::endl;
  }
}
/* Explicit specialisations for poq_distribution with double and float: */
template <class DerivedA, class DerivedB>
inline void quantile(const poq_distribution<double>& dist,
		     const Eigen::DenseBase<DerivedA>& x,
		     Eigen::DenseBase<DerivedB> const & y) {
  dist.quantile(x, y);
}
template <class DerivedA, class DerivedB>
inline void quantile(const poq_distribution<float>& dist,
		     const Eigen::DenseBase<DerivedA>& x,
		     Eigen::DenseBase<DerivedB> const & y) {
  dist.quantile(x, y);
}


template <class Dist, class DerivedA, class DerivedB>
void cdf(const Dist& dist,
	 const Eigen::DenseBase<DerivedA>& x,
	 Eigen::DenseBase<DerivedB> const & y_) {
  typename Eigen::DenseBase<DerivedB>& y =
    const_cast< typename Eigen::DenseBase<DerivedB>& >(y_);
  y.derived().resize(x.rows(), x.cols());
  for (int i=0; i < x.size(); ++i) {
    y.derived()[i] = cdf(dist, x.derived()[i]);
  }
}

template <class Dist, class DerivedA, class DerivedB>
void pdf(const Dist& dist,
	 const Eigen::DenseBase<DerivedA>& x,
	 Eigen::DenseBase<DerivedB> const & y_) {
  typename Eigen::DenseBase<DerivedB>& y =
    const_cast< typename Eigen::DenseBase<DerivedB>& >(y_);
  y.derived().resize(x.rows(), x.cols());
  for (int i=0; i < x.size(); ++i) {
    y.derived()[i] = pdf(dist, x.derived()[i]);
  }
}











//PoqEstimatorMethod

//template < int Method >

template<typename Scalar>
class poq_estimator {

  typedef typename poq_distribution<Scalar>::Vector Vector;
  typedef std::map<std::string, Scalar> SummaryMap;

  Vector x_;

  poq_distribution<Scalar> poq_;
  SummaryMap summary_;

  Scalar median() {
    return (x_[(int)((x_.size() - 1) / 2)] +
	    x_[x_.size() - 1 - (int)((x_.size() - 1) / 2)]) / 2.0;
  }
  Scalar mean_min() {
    return median() - x_.head(1 + (int)((x_.size() - 1) / 2)).mean();
  }
  Scalar mean_max() {
    return x_.tail(1 + (int)((x_.size() - 1) / 2)).mean() - median();
  }
  void summary_calc() {
    summary_.clear();
    summary_["median"] = median();
    summary_["m_min"] = mean_min();
    summary_["m_max"] = mean_max();
    summary_["m_sum"] = summary_["m_min"] + summary_["m_max"];
    summary_["m_rat"] = (summary_["m_max"] - summary_["m_min"]) / summary_["m_sum"];
    summary_["m_arat"] = std::fabs(summary_["m_rat"]);
  }
public:
  poq_estimator() {
    summary_.clear();
  }
  poq_estimator& data(const typename poq_distribution<Scalar>::Vector& x) {
    x_ = x;
    std::sort(x_.data(), x_.data()+x_.size());
    summary_calc();
    return *this;
  }
  const typename poq_distribution<Scalar>::Vector& data() {
    return x_;
  }
  void param(const typename poq_distribution<Scalar>::param_type& theta) {
    poq_.param(theta);
  }
  const typename poq_distribution<Scalar>::param_type& param() const {
    return poq_.param();
  }

  const SummaryMap* summary() {
    return &summary_;
  }
  /** 
   * Moment estimator
   * 
   * @param theta_ 
   * @param tayloroverride 
   */
  poq_estimator<Scalar>& moment(typename poq_distribution<Scalar>::Vector const & theta_,
				int tayloroverride = 0) {
    typename poq_distribution<Scalar>::Vector& theta =
      const_cast< typename poq_distribution<Scalar>::Vector& >(theta_);
    Scalar med = summary_["median"];
    Scalar m_sum = summary_["m_sum"];
    Scalar m_rat = summary_["m_rat"];
    Scalar m_arat = summary_["m_arat"];
    Scalar tmp = m_arat * std::log(2.0);
    if (((tmp < poq_.taylorLimitForward()) && (tayloroverride == 0))
	|| (tayloroverride < 0)) {
      tmp = 1.0 / std::log(2.0) / (1.0 - tmp / 2.0 * (1.0 - tmp / 3.0));
    } else {
      tmp = m_arat / (1.0 - std::exp(-tmp));
    }
    theta.resize(5);
    theta(0) = med;
    theta(1) = std::log(m_sum * (1.0 - m_arat) * tmp);
    theta(2) = logit(0.5 * (1.0 + m_rat / (tmp - 1.0)));
    theta(3) = m_arat;
    theta(4) = m_arat;
    return *this;
   };

  poq_estimator<Scalar>& moment(int tayloroverride = 0) {
    typename poq_distribution<Scalar>::Vector theta(5);
    moment(theta, tayloroverride);
    poq_.param(theta);
    return *this;
   };

  Scalar scoreLog(Vector const & score_) {
    return poq_.scoreLog(x_, score_);
  }
  Scalar scoreLog() {
    return poq_.scoreLog(x_);
  }
};

#endif
