#pragma once

#include <iostream>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <unsupported/Eigen/IterativeSolvers>

class MatrixReplacement;
using Eigen::SparseMatrix;

namespace Eigen {
    namespace internal {
        // MatrixReplacement looks-like a SparseMatrix, so let's inherits its traits:
        template<>
        struct traits<MatrixReplacement> : public Eigen::internal::traits<Eigen::SparseMatrix<double> >
        {};
    }
}

// Example of a matrix-free wrapper from a user type to Eigen's compatible type
// For the sake of simplicity, this example simply wrap a Eigen::SparseMatrix.
class MatrixReplacement : public Eigen::EigenBase<MatrixReplacement> {
public:
    // Required typedefs, constants, and method:
    typedef double Scalar;
    typedef double RealScalar;
    typedef int StorageIndex;
    enum {
        ColsAtCompileTime = Eigen::Dynamic,
        MaxColsAtCompileTime = Eigen::Dynamic,
        IsRowMajor = false
    };

    Index rows() const { return mp_mat->rows(); }
    Index cols() const { return mp_mat->cols(); }

    template<typename Rhs>
    Eigen::Product<MatrixReplacement, Rhs, Eigen::AliasFreeProduct> operator*(const Eigen::MatrixBase<Rhs>& x) const {
        return Eigen::Product<MatrixReplacement, Rhs, Eigen::AliasFreeProduct>(*this, x.derived());
    }

    // Custom API:
    MatrixReplacement() : mp_mat(0) {}

    void attachMyMatrix(const SparseMatrix<double>& mat) {
        mp_mat = &mat;
    }
    const SparseMatrix<double> my_matrix() const { return *mp_mat; }

private:
    const SparseMatrix<double>* mp_mat;
};


// Implementation of MatrixReplacement * Eigen::DenseVector though a specialization of internal::generic_product_impl:
namespace Eigen {
    namespace internal {

        template<typename Rhs>
        struct generic_product_impl<MatrixReplacement, Rhs, SparseShape, DenseShape, GemvProduct> // GEMV stands for matrix-vector
            : generic_product_impl_base<MatrixReplacement, Rhs, generic_product_impl<MatrixReplacement, Rhs> >
        {
            typedef typename Product<MatrixReplacement, Rhs>::Scalar Scalar;

            template<typename Dest>
            static void scaleAndAddTo(Dest& dst, const MatrixReplacement& lhs, const Rhs& rhs, const Scalar& alpha)
            {
                // This method should implement "dst += alpha * lhs * rhs" inplace,
                // however, for iterative solvers, alpha is always equal to 1, so let's not bother about it.
                assert(alpha == Scalar(1) && "scaling is not implemented");
                EIGEN_ONLY_USED_FOR_DEBUG(alpha);

                // Here we could simply call dst.noalias() += lhs.my_matrix() * rhs,
                // but let's do something fancier (and less efficient):
                for (Index i = 0; i < lhs.cols(); ++i)
                    dst += rhs(i) * lhs.my_matrix().col(i);
            }
        };

    }
}

namespace Eigen {
// This class simply warp a diagonal matrix as a Jacobi preconditioner.
// In the future such simple and generic wrapper should be shipped within Eigen itsel.
template <typename _Scalar>
class MyPreconditioner
{
    typedef _Scalar Scalar;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Vector;
    typedef typename Vector::Index Index;       // **
public:
    // this typedef is only to export the scalar type and compile-time dimensions to solve_retval
    typedef typename Vector::StorageIndex StorageIndex;
    enum {
        ColsAtCompileTime = Dynamic,
        MaxColsAtCompileTime = Dynamic
    };

    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> MatrixType;
    MyPreconditioner() : m_isInitialized(false) {}
    void setInvDiag(const Eigen::VectorXd& invdiag) {
        m_invdiag = invdiag;
        m_isInitialized = true;
    }
    Index rows() const { return m_invdiag.size(); }
    Index cols() const { return m_invdiag.size(); }

    template<typename MatType>
    MyPreconditioner& analyzePattern(const MatType&) { return *this; }

    template<typename MatType>
    MyPreconditioner& factorize(const MatType& mat) { return *this; }

    template<typename MatType>
    MyPreconditioner& compute(const MatType& mat) { return *this; }
    template<typename Rhs, typename Dest>
    void _solve_impl(const Rhs& b, Dest& x) const
    {
        x = m_invdiag.array() * b.array();
    }
    template<typename Rhs> inline const Solve<MyPreconditioner, Rhs>
        solve(const Eigen::MatrixBase<Rhs>& b) const
        {
            eigen_assert(m_isInitialized && "MyPreconditioner is not initialized.");
            eigen_assert(m_invdiag.size() == b.rows()
                && "MyPreconditioner::solve(): invalid number of rows of the right hand side matrix b");
            return Solve<MyPreconditioner, Rhs>(*this, b.derived());
        }

    ComputationInfo info() { return Success; }
protected:
    Vector m_invdiag;
    bool m_isInitialized;
};

template <typename _Scalar>
class MyPreconditioner2
{
    typedef _Scalar Scalar;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Vector;
    typedef typename Vector::Index Index;       // **
public:
    // this typedef is only to export the scalar type and compile-time dimensions to solve_retval
    typedef typename Vector::StorageIndex StorageIndex;
    enum {
        ColsAtCompileTime = Dynamic,
        MaxColsAtCompileTime = Dynamic
    };

    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> MatrixType;
    MyPreconditioner2() : m_isInitialized(false) {}
    void setInvDiag(const Eigen::VectorXd& approx) {
        m_approx = approx;
        m_isInitialized = true;
    }
    Index rows() const { return m_approx.size(); }
    Index cols() const { return m_approx.size(); }

    template<typename MatType>
    MyPreconditioner2& analyzePattern(const MatType&) { return *this; }

    template<typename MatType>
    MyPreconditioner2& factorize(const MatType& mat) { return *this; }

    template<typename MatType>
    MyPreconditioner2& compute(const MatType& mat) { return *this; }
    template<typename Rhs, typename Dest>
    void _solve_impl(const Rhs& b, Dest& x) const
    {
        //x = m_invdiag.array() * b.array();
        x = m_approx.array();
    }
    template<typename Rhs> inline const Solve<MyPreconditioner2, Rhs>
        solve(const Eigen::MatrixBase<Rhs>& b) const
        {
            eigen_assert(m_isInitialized && "MyPreconditioner2 is not initialized.");
            eigen_assert(m_approx.size() == b.rows()
                && "MyPreconditioner2::solve(): invalid number of rows of the right hand side matrix b");
            return Solve<MyPreconditioner2, Rhs>(*this, b.derived());
        }

        ComputationInfo info() { return Success; }
protected:
    Vector m_approx;
    bool m_isInitialized;
};

template <typename _Scalar>
class MyPreconditioner3
{
    typedef _Scalar Scalar;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Vector;
    typedef typename Vector::Index Index;       // **
public:
    // this typedef is only to export the scalar type and compile-time dimensions to solve_retval
    typedef typename Vector::StorageIndex StorageIndex;
    enum {
        ColsAtCompileTime = Dynamic,
        MaxColsAtCompileTime = Dynamic
    };

    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> MatrixType;
    MyPreconditioner3() : m_isInitialized(false) {}
    void setInvDiag(const Eigen::VectorXd& approx) {
        m_approx = approx;
        m_isInitialized = true;
    }
    Index rows() const { return m_approx.size(); }
    Index cols() const { return m_approx.size(); }

    template<typename MatType>
    MyPreconditioner3& analyzePattern(const MatType&) { return *this; }

    template<typename MatType>
    MyPreconditioner3& factorize(const MatType& mat) { return *this; }

    template<typename MatType>
    MyPreconditioner3& compute(const MatType& mat) { return *this; }
    template<typename Rhs, typename Dest>
    void _solve_impl(const Rhs& b, Dest& x) const
    {
        //x = m_invdiag.array() * b.array();
        x = b.array();
    }
    template<typename Rhs> inline const Solve<MyPreconditioner3, Rhs>
        solve(const Eigen::MatrixBase<Rhs>& b) const
        {
            eigen_assert(m_isInitialized && "MyPreconditioner3 is not initialized.");
            eigen_assert(m_approx.size() == b.rows()
                && "MyPreconditioner3::solve(): invalid number of rows of the right hand side matrix b");
            return Solve<MyPreconditioner3, Rhs>(*this, b.derived());
        }

        ComputationInfo info() { return Success; }
protected:
    Vector m_approx;
    bool m_isInitialized;
};

template <typename _Scalar>
class ILU_0_Preconditioner
{
    typedef _Scalar Scalar;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Vector;
    typedef typename Vector::Index Index;       // **
public:
    // this typedef is only to export the scalar type and compile-time dimensions to solve_retval
    typedef typename Vector::StorageIndex StorageIndex;
    enum {
        ColsAtCompileTime = Dynamic,
        MaxColsAtCompileTime = Dynamic
    };

    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> MatrixType;
    ILU_0_Preconditioner() : m_isInitialized(false) {}
    void setLUmatrix(const Eigen::SparseMatrix< double >& LU, int eq_num_) {
        eq_num = eq_num_;
        m_LU = LU;
        m_isInitialized = true;
    }
    Index rows() const { return m_LU.rows(); }
    Index cols() const { return m_LU.cols(); }

    template<typename MatType>
    ILU_0_Preconditioner& analyzePattern(const MatType&) { return *this; }

    template<typename MatType>
    ILU_0_Preconditioner& factorize(const MatType& mat) { return *this; }

    template<typename MatType>
    ILU_0_Preconditioner& compute(const MatType& mat) { return *this; }
    template<typename Rhs, typename Dest>
    void _solve_impl(const Rhs& b, Dest& x) const
    {
        Vector y(rows()), x_(rows());
        int j_start, j_finish;

        y(0) = b(0);
        for (int i = 1; i < rows(); ++i)
        {
            j_start = max(0, i / eq_num - 1) * eq_num;
            y(i) = b(i);
            for (int j = j_start; j < i; ++j) {
                y(i) -= y(j) * m_LU.coeff(i, j);
            }
        }
        x_(rows() - 1) = y(rows() - 1) / m_LU.coeff(rows() - 1, rows() - 1);
        for (int i = rows() - 1 - 1; i >= 0 ; --i)
        {
            j_finish = (i / eq_num - 1 + 3) * eq_num;
            if (j_finish > rows()) j_finish = rows();

            x_(i) = y(i);
            for (int j = i + 1; j < j_finish; ++j) {
                x_(i) -= x_(j) * m_LU.coeff(i, j);
            }
            x_(i) /= m_LU.coeff(i, i);
        }
        x = x_.array();
    }
    template<typename Rhs> inline const Solve<ILU_0_Preconditioner, Rhs>
        solve(const Eigen::MatrixBase<Rhs>& b) const
        {
            eigen_assert(m_isInitialized && "ILU_0_Preconditioner is not initialized.");
            eigen_assert(m_LU.rows() == b.rows()
                && "ILU_0_Preconditioner::solve(): invalid number of rows of the right hand side matrix b");
            return Solve<ILU_0_Preconditioner, Rhs>(*this, b.derived());
        }

        ComputationInfo info() { return Success; }
protected:
    SparseMatrix< double > m_LU;
    int eq_num;
    bool m_isInitialized;
};

}

//int main()
//{
//    int n = 10;
//    Eigen::SparseMatrix<double> S = Eigen::MatrixXd::Random(n, n).sparseView(0.5, 1);
//    S = S.transpose() * S;
//
//    MatrixReplacement A;
//    A.attachMyMatrix(S);
//
//    Eigen::VectorXd b(n), x;
//    b.setRandom();
//
//    // Solve Ax = b using various iterative solver with matrix-free version:
//    {
//        Eigen::ConjugateGradient<MatrixReplacement, Eigen::Lower | Eigen::Upper, Eigen::IdentityPreconditioner> cg;
//        cg.compute(A);
//        x = cg.solve(b);
//        std::cout << "CG:       #iterations: " << cg.iterations() << ", estimated error: " << cg.error() << std::endl;
//    }
//
//    {
//        Eigen::BiCGSTAB<MatrixReplacement, Eigen::IdentityPreconditioner> bicg;
//        bicg.compute(A);
//        x = bicg.solve(b);
//        std::cout << "BiCGSTAB: #iterations: " << bicg.iterations() << ", estimated error: " << bicg.error() << std::endl;
//    }
//
//    {
//        Eigen::GMRES<MatrixReplacement, Eigen::IdentityPreconditioner> gmres;
//        gmres.compute(A);
//        x = gmres.solve(b);
//        std::cout << "GMRES:    #iterations: " << gmres.iterations() << ", estimated error: " << gmres.error() << std::endl;
//    }
//
//    {
//        Eigen::DGMRES<MatrixReplacement, Eigen::IdentityPreconditioner> gmres;
//        gmres.compute(A);
//        x = gmres.solve(b);
//        std::cout << "DGMRES:   #iterations: " << gmres.iterations() << ", estimated error: " << gmres.error() << std::endl;
//    }
//
//    {
//        Eigen::MINRES<MatrixReplacement, Eigen::Lower | Eigen::Upper, Eigen::IdentityPreconditioner> minres;
//        minres.compute(A);
//        x = minres.solve(b);
//        std::cout << "MINRES:   #iterations: " << minres.iterations() << ", estimated error: " << minres.error() << std::endl;
//    }
//}