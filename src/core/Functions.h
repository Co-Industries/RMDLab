#ifndef RMD_Functions_h
#define RMD_Functions_h

#include <Corrade/Containers/Array.h>
#include <Magnum/Magnum.h>
#include <Magnum/Math/Color.h>

namespace Magnum
{
    using namespace Math::Literals;

    // ? QEq qeq.F90
    void QEq(); /* Two vector electronegativity equilization routine */
    Double dot_product(const Int &mode);

    // ? FORCE pot.F90
    void FORCE();
    void ENbond(); /* Calculates the energy and the forces due to the Van der Waals and Coulomb terms */
    void Ebond();
    void ForceBbo(const std::size_t &i, const std::size_t &j, const std::size_t &b, const Vector3d &coeff);
    void Elnpr();

    /*
    *To find out hydrogen bonding combinations, <vnhbp> (one atom parameter, 2nd row, 8th column) is used to 
    *identify whether an atoms is hydrogen or not, <vnhbp>=1 for H, <vnhbp>=2 for O,N,S and <vnhbp>=0 for others. 
    *If a 3atom combination doesn't consists of H,O,N,S, the hydrogen bonding doens't exist.
    *Among the three atoms, the center atom is always a hydrogen. In the original code, atom IDs are given 
    *in such as j2 <-> j1(H) <=> i2. <-> is bonding interaction, <=> is non-bonding interaction.
    *If the magnitude of the donor/acceptor bond is less than 1.d-2, the hydrogen bonding doesn't exit.
    *If the interatomic distance between j2 <-> i2 is larger than 10[A], the interaction doesn't exit. 
    */
    void Ehb();

    /*
    *Derivative of BOij using new bond order definition. Only difference is that sigma BO
    *prime is replaced with full BOp. The derivative of BO becomes a bit simpler due to
    *the new definition.
    */
    void ForceB(const std::size_t &i, const std::size_t &j, const std::size_t &b, const Double &coeff);

    // *derivative of <cos_ijk>
    void ForceA3(const Double &coeff, const std::size_t &i, const std::size_t &j, const std::size_t &k, const Vector3d &da0, const Vector3d &da1, const Double &da0_0, const Double &da1_0);

    // ? bo.F90
    void BOCALC();
    void BOPRIM(); /* Calculates the BOp(0:3,i,j) and the deltap(i) */
    void BOFULL(); /* Calculates the Bond Order and its derivatives */
}

#endif