#ifndef MASS_SPRING_HPP
#define MASS_SPRING_HPP

#include <nonlinfunc.hpp>
#include <timestepper.hpp>
#include <autodiff_dynamic.hpp>

using namespace ASC_ode;

#include <vector.hpp>
using namespace nanoblas;


template <int D>
class Mass
{
public:
  double mass;
  Vec<D> pos;
  Vec<D> vel = 0.0;
  Vec<D> acc = 0.0;
};


template <int D>
class Fix
{
public:
  Vec<D> pos;
};


class Connector
{
public:
  enum CONTYPE { FIX=1, MASS=2 };
  CONTYPE type;
  size_t nr;
};


std::ostream & operator<< (std::ostream & ost, const Connector & con)
{
  ost << "type = " << int(con.type) << ", nr = " << con.nr;
  return ost;
}

class Spring
{
public:
  double length;  
  double stiffness;
  std::array<Connector,2> connectors;
};

class Joint 
{
public:
  double length;
  std::array<Connector,2> connectors;
};

template <int D>
class MassSpringSystem
{
  std::vector<Fix<D>> m_fixes;
  std::vector<Mass<D>> m_masses;
  std::vector<Spring> m_springs;
  std::vector<Joint> m_joints;
  Vec<D> m_gravity=0.0;
public:
  void setGravity (Vec<D> gravity) { m_gravity = gravity; }
  Vec<D> getGravity() const { return m_gravity; }

  Connector addFix (Fix<D> p)
  {
    m_fixes.push_back(p);
    return { Connector::FIX, m_fixes.size()-1 };
  }

  Connector addMass (Mass<D> m)
  {
    m_masses.push_back (m);
    return { Connector::MASS, m_masses.size()-1 };
  }
  
  size_t addSpring (Spring s) 
  {
    m_springs.push_back (s); 
    return m_springs.size()-1;
  }
  size_t addJoint (Joint j) 
  {
    m_joints.push_back (j); 
    return m_joints.size()-1;
  }

  auto & fixes() { return m_fixes; } 
  auto & masses() { return m_masses; } 
  auto & springs() { return m_springs; }
  auto & joints() { return m_joints; }

  void getState (VectorView<> values, VectorView<> dvalues, VectorView<> ddvalues)
  {
    auto valmat = values.asMatrix(m_masses.size(), D);
    auto dvalmat = dvalues.asMatrix(m_masses.size(), D);
    auto ddvalmat = ddvalues.asMatrix(m_masses.size(), D);

    for (size_t i = 0; i < m_masses.size(); i++)
      {
        valmat.row(i) = m_masses[i].pos;
        dvalmat.row(i) = m_masses[i].vel;
        ddvalmat.row(i) = m_masses[i].acc;
      }
  }

  void setState (VectorView<> values, VectorView<> dvalues, VectorView<> ddvalues)
  {
    auto valmat = values.asMatrix(m_masses.size(), D);
    auto dvalmat = dvalues.asMatrix(m_masses.size(), D);
    auto ddvalmat = ddvalues.asMatrix(m_masses.size(), D);

    for (size_t i = 0; i < m_masses.size(); i++)
      {
        m_masses[i].pos = valmat.row(i);
        m_masses[i].vel = dvalmat.row(i);
        m_masses[i].acc = ddvalmat.row(i);
      }
  }
};

template <int D>
std::ostream & operator<< (std::ostream & ost, MassSpringSystem<D> & mss)
{
  ost << "fixes:" << std::endl;
  for (auto f : mss.fixes())
    ost << f.pos << std::endl;

  ost << "masses: " << std::endl;
  for (auto m : mss.masses())
    ost << "m = " << m.mass << ", pos = " << m.pos << std::endl;

  ost << "springs: " << std::endl;
  for (auto sp : mss.springs())
    ost << "length = " << sp.length << ", stiffness = " << sp.stiffness
        << ", C1 = " << sp.connectors[0] << ", C2 = " << sp.connectors[1] << std::endl;


  ost << "joints: " << std::endl;
  for (auto j : mss.joints())
    ost << "length = " << j.length
        << ", C1 = " << j.connectors[0] << ", C2 = " << j.connectors[1] << std::endl;

  return ost;
}


template <int D>
class MSS_Function : public NonlinearFunction
{
  MassSpringSystem<D> & mss;
public:
  MSS_Function (MassSpringSystem<D> & _mss)
    : mss(_mss) { }

  virtual size_t dimX() const override { return D*mss.masses().size() + mss.joints().size(); }
  virtual size_t dimF() const override{ return D*mss.masses().size() +  mss.joints().size(); }

  virtual void evaluate(VectorView<double> x, VectorView<double> f) const override {
    evaluateGeneric(x, f);
  }

  template<typename T>
  void evaluateGeneric(VectorView<T> x, VectorView<T> f) const {
    f = 0.0;

    int lamdacounter = mss.joints().size();

    auto xmat = x.asMatrix(mss.masses().size(), D);
    auto xlambda = x.range(D*mss.masses().size(), lamdacounter);
    auto fmat = f.asMatrix(mss.masses().size(), D);
    auto flambda = f.range(D*mss.masses().size(), lamdacounter);

    // gravity force
    for (size_t i = 0; i < mss.masses().size(); i++)
      fmat.row(i) = mss.masses()[i].mass*mss.getGravity();

    // spring forces
    for (auto spring : mss.springs())
    {
        auto [c1, c2] = spring.connectors;
        Vec<D, T> p1, p2;
        if (c1.type == Connector::FIX)
          p1 = mss.fixes()[c1.nr].pos; // implicit cast double->T
        else
          p1 = xmat.row(c1.nr);
        if (c2.type == Connector::FIX)
          p2 = mss.fixes()[c2.nr].pos;
        else
          p2 = xmat.row(c2.nr);

        Vec<D, T> p2MinusP1 = p2-p1;
        T dist = norm(p2MinusP1);

        T force = spring.stiffness * (dist - spring.length);
        Vec<D, T> dir12 = (1.0/dist) * p2MinusP1;

        if (c1.type == Connector::MASS)
          fmat.row(c1.nr) += force*dir12;
        if (c2.type == Connector::MASS)
          fmat.row(c2.nr) -= force*dir12;
      }

      // joint part
      for (size_t i=0; i<lamdacounter; i ++)
      {
        Joint joint = mss.joints()[i];
        auto [c1, c2] = joint.connectors;
        Vec<D, T> p1, p2;
        if (c1.type == Connector::FIX && c2.type == Connector::FIX)
          throw std::invalid_argument("Both connectors of a joint cannot be fixed.");

        if (c1.type == Connector::FIX)
          p1 = mss.fixes()[c1.nr].pos;
        else
          p1 = xmat.row(c1.nr);
        if (c2.type == Connector::FIX)
          p2 = mss.fixes()[c2.nr].pos;
        else
          p2 = xmat.row(c2.nr);

        Vec<D, T> p1MinusP2 = p1-p2;
        Vec<D, T> force_vec;
        for (int j = 0; j < D; j++) {
          force_vec(j) = p1MinusP2(j) * 2.0 * xlambda(i);
        }

        if (c1.type == Connector::MASS)
          fmat.row(c1.nr) += force_vec;
        if (c2.type == Connector::MASS)
          fmat.row(c2.nr) -= force_vec;

        T temp = 0.0;
        for (size_t k = 0; k < D; k++)
          temp += (p1(k) - p2(k)) * (p1(k) - p2(k));

        flambda(i) = temp - joint.length * joint.length;
      }

      for (size_t i = 0; i < mss.masses().size(); i++) {
        for (size_t j = 0; j < D; j++) {
          fmat(i, j) = fmat(i, j) / mss.masses()[i].mass;
        }
      }
    }

    virtual void evaluateDeriv(VectorView<double> x, MatrixView<double> df) const override {
      Vector<AutoDiffDynamic<>> xs(x.size());
      Vector<AutoDiffDynamic<>> fs(x.size());
      for (int i = 0; i < x.size(); i++) {
        xs(i) = AutoDiffDynamic(x(i), x.size()); // <--- FIX: Assign to index 'i'
        xs(i).deriv()[i] = 1.0;
      }
      evaluateGeneric(xs, fs);
      for (size_t r = 0; r < dimF(); r++) {
        for (size_t c = 0; c < dimX(); c++) {
          df(r, c) = fs[r].deriv()[c];
        }
      }
    }

};

#endif