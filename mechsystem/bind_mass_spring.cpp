#include <sstream>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include "mass_spring.hpp"
#include "Newmark.hpp"

namespace py = pybind11;

PYBIND11_MAKE_OPAQUE(std::vector<Mass<3>>);
PYBIND11_MAKE_OPAQUE(std::vector<Fix<3>>);
PYBIND11_MAKE_OPAQUE(std::vector<Spring>);

PYBIND11_MODULE(mass_spring, m) {
    m.doc() = "mass-spring-system simulator"; 

    py::class_<Mass<2>> (m, "Mass2d")
          .def_property("mass",
                    [](Mass<2> & m) { return m.mass; },
                    [](Mass<2> & m, double mass) { m.mass = mass; })
      .def_property_readonly("pos",
                             [](Mass<2> & m) { return m.pos.data(); });

      ;
      
    m.def("Mass", [](double m, std::array<double,2> p)
    {
      return Mass<2>{m, { p[0], p[1] }};
    });

    
    py::class_<Mass<3>> (m, "Mass3d")
      .def_property("mass",
                    [](Mass<3> & m) { return m.mass; },
                    [](Mass<3> & m, double mass) { m.mass = mass; })
      .def_property_readonly("pos",
                             [](Mass<3> & m) { return m.pos.data(); });
    ;

    m.def("Mass", [](double m, std::array<double,3> p, std::array<double,3> v)
    {
      return Mass<3>{m, { p[0], p[1], p[2] }, { v[0], v[1], v[2] }};
    });

    m.def("Mass", [](double m, std::array<double,3> p, std::array<double,3> v)
    {
      return Mass<3>{m, { p[0], p[1], p[2] }, { v[0], v[1], v[2] }};
    });


    py::class_<Fix<2>> (m, "Fix2d")
      .def_property_readonly("pos",
                             [](Fix<2> & f) { return f.pos.data(); });

    m.def("Fix", [](std::array<double,2> p)
    {
      return Fix<2>{ { p[0], p[1] } };
    });


    
    py::class_<Fix<3>> (m, "Fix3d")
      .def_property_readonly("pos",
                             [](Fix<3> & f) { return f.pos.data(); });
    
    m.def("Fix", [](std::array<double,3> p)
    {
      return Fix<3>{ { p[0], p[1], p[2] } };
    });

    py::class_<Connector> (m, "Connector");

    py::class_<Spring> (m, "Spring")
      .def(py::init<double, double, std::array<Connector,2>>())
      .def_property_readonly("connectors",
                             [](Spring & s) { return s.connectors; })
      ;

    py::class_<Joint> (m, "Joint")
        .def(py::init<double, std::array<Connector,2>>())
        .def_property_readonly("connectors",
                               [](Joint & j) { return j.connectors; })
        ;

    
    py::bind_vector<std::vector<Mass<3>>>(m, "Masses3d");
    py::bind_vector<std::vector<Fix<3>>>(m, "Fixes3d");
    py::bind_vector<std::vector<Spring>>(m, "Springs");        
    
    
    py::class_<MassSpringSystem<2>> (m, "MassSpringSystem2d")
      .def(py::init<>())
      .def("add", [](MassSpringSystem<2> & mss, Mass<2> m) { return mss.addMass(m); })
      ;
      
        
    py::class_<MassSpringSystem<3>> (m, "MassSpringSystem3d")
      .def(py::init<>())
      .def("__str__", [](MassSpringSystem<3> & mss) {
        std::stringstream sstr;
        sstr << mss;
        return sstr.str();
      })
      .def_property("gravity", [](MassSpringSystem<3> & mss) { return mss.getGravity(); },
                    [](MassSpringSystem<3> & mss, std::array<double,3> g) { mss.setGravity(Vec<3>{g[0],g[1],g[2]}); })
      .def("add", [](MassSpringSystem<3> & mss, Mass<3> m) { return mss.addMass(m); })
      .def("add", [](MassSpringSystem<3> & mss, Fix<3> f) { return mss.addFix(f); })
      .def("add", [](MassSpringSystem<3> & mss, Spring s) { return mss.addSpring(s); })
      .def("add", [](MassSpringSystem<3> & mss, Joint j) { return mss.addJoint(j); })
      .def_property_readonly("masses", [](MassSpringSystem<3> & mss) -> auto& { return mss.masses(); })
      .def_property_readonly("fixes", [](MassSpringSystem<3> & mss) -> auto& { return mss.fixes(); })
      .def_property_readonly("springs", [](MassSpringSystem<3> & mss) -> auto& { return mss.springs(); })
      .def_property_readonly("joints", [](MassSpringSystem<3> & mss) -> auto& { return mss.joints(); })
      .def("__getitem__", [](MassSpringSystem<3> mss, Connector & c) {
        if (c.type==Connector::FIX) return py::cast(mss.fixes()[c.nr]);
        else return py::cast(mss.masses()[c.nr]);
      })
      
      .def("getState", [] (MassSpringSystem<3> & mss) {
        Vector<> x(3*mss.masses().size());
        Vector<> dx(3*mss.masses().size());
        Vector<> ddx(3*mss.masses().size());
        mss.getState (x, dx, ddx);
        return std::vector<double>(x);
      })

      // call solver
      .def("simulate", [](MassSpringSystem<3> & mss, double tend, size_t steps) {
        Vector<> x(3*mss.masses().size()+ mss.joints().size());
        Vector<> dx(3*mss.masses().size()+ mss.joints().size());
        Vector<> ddx(3*mss.masses().size()+ mss.joints().size());
        mss.getState (x, dx, ddx);

        Vector<double> temp (mss.masses().size()*3 + mss.joints().size());
        temp = 1e-12; // Small regularization for Lagrange multipliers (lower = more accurate constraints)
        for (size_t i = 0; i < mss.masses().size(); i++){
            temp(3*i) = mss.masses()[i].mass;
            temp(3*i+1) = mss.masses()[i].mass;
            temp(3*i+2) = mss.masses()[i].mass;
        }

        auto mss_func = std::make_shared<MSS_Function<3>> (mss);
        auto mass = std::make_shared<DiagMatrixFunction> (temp);

        SolveODE_Alpha(tend, steps, 0.8, x, dx, ddx, mss_func, mass);

        mss.setState (x, dx, ddx);  
    });
    // .def("updateLengthConstraints", [](MassSpringSystem<3> & mss, Vector<> & f, const Vector<> & x) {
    //     const size_t D = 3;
    //     // length constraints
    //     for (size_t i = 0; i < mss.joints().size(); i++){
    //         for (size_t k = 0; k < D; k++)
    //           {
    //             auto idx_0 = mss.joints()[i].connectors[0].nr;
    //             if (mss.joints()[i].connectors[0].type == Connector::FIX)
    //               f(D*mss.masses().size()+i) += std::pow(x(i*D+k) - mss.fixes()[idx_0].pos(k), 2);
    //             else
    //               f(D*mss.masses().size()+i) += std::pow(x(i*D+k) - mss.masses()[idx_0].pos(k), 2);
    //           }
    //         f(D*mss.masses().size()+i) -= std::pow(mss.joints()[i].length, 2);
    //     }
    //     // d/dx <lambda, g(x)>
    //     for (size_t i = 0; i <mss.joints().size(); i++)
    //        for (size_t j =0; j < D; j++)
    //           f(D*i+j) += 2*(mss.masses()[i].pos(j) - x(D*i+j))*x(D*mss.masses().size()+i);
        
    //     // d/dx <lambda, g(x)> = lambda * gradient of constraint
    //     // For constraint g = |p1-p2|^2 - L^2, gradient is 2*(p1-p2)
    //     for (size_t i = 0; i < mss.joints().size(); i++)
    //     {
    //         double lambda = x(D*mss.masses().size()+i);  // Extract lambda from x
    //         auto [c1, c2] = mss.joints()[i].connectors;
            
    //         // Extract positions from x vector (not stored positions!)
    //         Vec<D> p1, p2;
    //         if (c1.type == Connector::FIX)
    //             p1 = mss.fixes()[c1.nr].pos;
    //         else
    //             for (size_t j = 0; j < D; j++)
    //                 p1(j) = x(D*c1.nr+j);  // From x vector at mass c1.nr
                    
    //         if (c2.type == Connector::FIX)
    //             p2 = mss.fixes()[c2.nr].pos;
    //         else
    //             for (size_t j = 0; j < D; j++)
    //                 p2(j) = x(D*c2.nr+j);  // From x vector at mass c2.nr
        
    //         // Apply constraint forces: f += lambda * gradient_of_g
    //         // gradient_of_g w.r.t p1 is 2*(p1-p2), w.r.t p2 is 2*(p2-p1)
    //         if (c1.type == Connector::MASS)
    //             for (size_t j = 0; j < D; j++)
    //                 f(D*c1.nr+j) += 2*lambda*(p1(j) - p2(j));
                    
    //         if (c2.type == Connector::MASS)
    //             for (size_t j = 0; j < D; j++)
    //                 f(D*c2.nr+j) += 2*lambda*(p2(j) - p1(j));
    //     }

    //     // Constraint equations: g(x) = |p1-p2|^2 - L^2 = 0
    //     for (size_t i = 0; i < mss.joints().size(); i++)
    //     {
    //         auto [c1, c2] = mss.joints()[i].connectors;
            
    //         // Extract BOTH connector positions from x vector
    //         Vec<D> p1, p2;
    //         if (c1.type == Connector::FIX)
    //             p1 = mss.fixes()[c1.nr].pos;
    //         else
    //             for (size_t j = 0; j < D; j++)
    //                 p1(j) = x(D*c1.nr+j);  // Mass position from x
                    
    //         if (c2.type == Connector::FIX)
    //             p2 = mss.fixes()[c2.nr].pos;
    //         else
    //             for (size_t j = 0; j < D; j++)
    //                 p2(j) = x(D*c2.nr+j);  // Mass position from x
        
    //         // Compute constraint: distance^2 - length^2
    //         double dist_sq = 0;
    //         for (size_t k = 0; k < D; k++)
    //             dist_sq += std::pow(p1(k) - p2(k), 2);
                
    //         f(D*mss.masses().size()+i) = dist_sq - std::pow(mss.joints()[i].length, 2);
    //     }
    // });
}
