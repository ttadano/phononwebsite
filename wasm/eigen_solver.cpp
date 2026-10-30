#include <complex>

#include <Eigen/Core>
#include <Eigen/Eigenvalues>

extern "C" int solve_hermitian_eigen(
    int size,
    const double* matrix_interleaved,
    double* eigenvalues_out,
    double* eigenvectors_out
) {
    using Complex = std::complex<double>;
    using MatrixXcd = Eigen::Matrix<Complex, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;
    using VectorXd = Eigen::VectorXd;

    if (size <= 0 || !matrix_interleaved || !eigenvalues_out || !eigenvectors_out) {
        return 2;
    }

    Eigen::Map<const MatrixXcd> matrix(
        reinterpret_cast<const Complex*>(matrix_interleaved),
        size,
        size
    );
    MatrixXcd hermitian = matrix;

    Eigen::SelfAdjointEigenSolver<MatrixXcd> solver(hermitian, Eigen::ComputeEigenvectors);
    if (solver.info() != Eigen::Success) {
        return 1;
    }

    Eigen::Map<VectorXd>(eigenvalues_out, size) = solver.eigenvalues();
    Eigen::Map<MatrixXcd>(reinterpret_cast<Complex*>(eigenvectors_out), size, size) = solver.eigenvectors();
    return 0;
}
