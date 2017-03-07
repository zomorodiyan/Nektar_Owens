///////////////////////////////////////////////////////////////////////////////
//
// File ViscoElasticConstitutive.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// License for the specific language governing rights and limitations under
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
//
// Description: Basic Class for Calculation ViscoElastic Stress,
//  OldRoyd-B, FENE-P, FENE-CR, Giesekus models
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_ViscoElasticConstitutive_H
#define NEKTAR_SOLVERS_ViscoElasticConstitutive_H

#include <HomogeneousBloodFlowSolver/EquationSystems/IncNavierStokes.h>
#include <MultiRegions/GlobalLinSysDirectStaticCond.h>

#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LibUtilities/FFT/NektarFFT.h>  // for NektarFFTSharedPtr
#include <SpatialDomains/MeshGraph.h>   // for MeshGraphSharedPtr
#include <MultiRegions/ExpList.h>       // for ExpListSharedPtr

///


namespace Nektar
{     

    /*   enum ViscoElasticType */
    /*  { */
    /*    eNoViscoElasticType, */
    /*    eOldroydB, */
    /*    eHomogeneousBloodModel, */
    /*    eUCM, */
    /*    eEquationViscoElasticTypeSize, */
    /*  }; */


    /*  // Keep this consistent with the enums in EquationType. */
    /*   //   const std::string kConstitutiveModelStr[] = */
    /* const std::string kEquationViscoElasticTypeStr [] =  */
    /*  { */
    /*      "NoType", */
    /*      "OldRoydB", */
    /* 	"HomogeneousBloodModel", */
    /*      "UCM" */
    /*  }; */

    enum   ViscoElasticTreatmentType
    {
        eNOViscoElasticTreatmentType,
        eImplicitExplicit,
        eFullExplicit,
        eViscoElasticTreatmentTypeSize,


    };

    // Keep this consistent with the enums in ViscoElasticTreatmentType
    const std::string kViscoElasticTreatmentTypeStr [] = 
    {
        "NOType",
        "ImplicitExplicit",
        "FullExplicit"

    };


    /* static NekDouble kHighOrderBCsExtrapolation[][3] = {{ 1.0,  0.0, 0.0}, */
    /*                                                      { 2.0, -1.0, 0.0}, */
    /*                                                      { 3.0, -3.0, 1.0}}; */


    class ViscoElasticConstitutive: public IncNavierStokes
    {


        public:

            /// Creates an instance of this class
            static SolverUtils::EquationSystemSharedPtr create(
                    const LibUtilities::SessionReaderSharedPtr& pSession) {
                SolverUtils::EquationSystemSharedPtr p = MemoryManager<ViscoElasticConstitutive>::AllocateSharedPtr(pSession);
                p->InitObject();
                return p;
            }

            /// Name of class
            static std::string className;


            //j ViscoElasticType m_viscoelasticType;   ///ViscoElastic Type
            ViscoElasticTreatmentType m_viscoelastictreatmentType; // ViscoElasticTreatment 

            // Destructor 

            virtual ~ViscoElasticConstitutive();

            //Constructor 

            ViscoElasticConstitutive(const LibUtilities::SessionReaderSharedPtr& pSession);




            void v_InitObject();

            void EvaluateViscoElasticStress_OldroydB_FullExpTreatment_RHS(const Array<OneD, const Array<OneD, NekDouble> > &inarray, 
                    Array<OneD, Array<OneD, NekDouble> > &outarray, 
                    const NekDouble time);

//-------------------------------------------Me.---------------------------------------------
            void EvaluateAggregateSize_FullExpTreatment_RHS(const Array<OneD, const Array<OneD, NekDouble> > &inarray,
                    Array<OneD, Array<OneD, NekDouble> > &outarray,
                    const NekDouble time);
//-------------------------------------------Me.---------------------------------------------
            void EvaluateViscoElasticStress_OldroydB_FullExpTreatment(const Array<OneD, const Array<OneD, NekDouble> > &inarray, 
                    Array<OneD, Array<OneD, NekDouble> > &outarray, 
                    const NekDouble time,
                    const NekDouble aii_Dt);

//-------------------------------------------Me.---------------------------------------------
            void EvaluateAggregateSize_FullExpTreatment(const Array<OneD, const Array<OneD, NekDouble> > &inarray, 
                    Array<OneD, Array<OneD, NekDouble> > &outarray, 
                    const NekDouble time,
                    const NekDouble aii_Dt);

//-------------------------------------------Me.---------------------------------------------
            void EvaluateViscoElasticStress_OldroydB_Homo_FullExpTreatment_RHS(const Array<OneD, const Array<OneD, NekDouble> > &inarray,
                    Array<OneD, Array<OneD, NekDouble> > &outarray,
                    const NekDouble time);

            void EvaluateViscoElasticStress_OldroydB__Homo_FullExpTreatment(const Array<OneD, const Array<OneD, NekDouble> > &inarray,
                    Array<OneD, Array<OneD, NekDouble> > &outarray,
                    const NekDouble time,
                    const NekDouble aii_Dt);


            void SetUpViscoElasticForcing(Array<OneD, Array<OneD, NekDouble> > &Forcing);

            void SetUpViscoElasticForcing_Weak(Array<OneD, Array<OneD, NekDouble> > &Forcing);


            void AddStressTimesNormalToVelocityNeumannBC(const Array<OneD, const Array<OneD, NekDouble> > &stressin,
                    Array<OneD, Array<OneD, NekDouble> > &Weakoutarray);


            void EvaluateViscoElasticStress_OldroydB_ImpExpTreatment_ExplicitPart(const Array<OneD, const Array<OneD, NekDouble> > &inarray, 
                    Array<OneD, Array<OneD, NekDouble> > &outarray, 
                    const NekDouble time);

            void EvaluateViscoElasticStress_OldroydB_ImpExpTreatment_Implicitpart(const Array<OneD, const Array<OneD, NekDouble> > &inarray,
                    Array<OneD, Array<OneD, NekDouble> > &outarrayCoeffs,
                    const NekDouble time,
                    const NekDouble aii_Dt);

            void EvaluateViscoElasticStress_OldroydB_ImpExpTreatment_ExplicitPart_Weak(const Array<OneD, const Array<OneD, NekDouble> > &inarray, 
                    Array<OneD, Array<OneD, NekDouble> > &outarray, 
                    const NekDouble time);
    };

}// end of namespace
#endif //NEKTAR_SOLVERS_ViscoElasticConstitutive_H
