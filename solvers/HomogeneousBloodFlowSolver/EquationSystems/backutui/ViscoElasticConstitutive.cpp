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
//  OldRoyd-B,  model
//
///////////////////////////////////////////////////////////////////////////////

#include <ViscoElasticFlowSolver/EquationSystems/ViscoElasticConstitutive.h>
#include <LibUtilities/BasicUtils/Timer.h>
#include <MultiRegions/ContField3DHomogeneous1D.h>
#include <MultiRegions/ContField3DHomogeneous2D.h>
//#include <Auxiliary/EquationSystem.h>
//**************
#include <MultiRegions/GlobalLinSysDirectStaticCond.h>

namespace Nektar
{

    string ViscoElasticConstitutive::className = SolverUtils::GetEquationSystemFactory().RegisterCreatorFunction(" ViscoElasticConstitutive", ViscoElasticConstitutive::create);

    ViscoElasticConstitutive::ViscoElasticConstitutive(const LibUtilities::SessionReaderSharedPtr& pSession):
        IncNavierStokes(pSession)
    {
    }

    void ViscoElasticConstitutive::v_InitObject()
    {

        int i;
        int physTot=m_fields[0]->GetTotPoints();
        IncNavierStokes::v_InitObject();
        m_CoeffState = MultiRegions::eLocal;

        if(m_equationType == eUnsteadyViscoElastic)
        {
            for(i = 0; i < (int) eViscoElasticTreatmentTypeSize; ++i)
            {
                bool match;
                m_session->MatchSolverInfo("VISCOELASTICTREATMENT",kViscoElasticTreatmentTypeStr[i],match,false); cout<< "test"<<kViscoElasticTreatmentTypeStr[i]<< endl; if(match)
                {
                    m_viscoelastictreatmentType = (ViscoElasticTreatmentType)i; 
                    break;
                }
            }

            ASSERTL0(i !=eViscoElasticTreatmentTypeSize ,"VISCOELASTICTREATMENT not found in SOLVERINFO section");
            //----------------------------------------------Me.------------------------------------------------
            switch(m_viscoelasticType)
            {
                case eHomogeneousBloodModel:

                    m_integrationOps_aggregatesize.DefineOdeRhs(&ViscoElasticConstitutive::EvaluateAggregateSize_FullExpTreatment_RHS, this);

                    m_integrationOps_aggregatesize.DefineImplicitSolve(&ViscoElasticConstitutive::EvaluateAggregateSize_FullExpTreatment, this);

                    m_integrationOps_viscoelasticstress.DefineOdeRhs(&ViscoElasticConstitutive::EvaluateViscoElasticStress_OldroydB_FullExpTreatment_RHS, this);

                    m_integrationOps_viscoelasticstress.DefineImplicitSolve(&ViscoElasticConstitutive::EvaluateViscoElasticStress_OldroydB_FullExpTreatment, this);


                    break;
                    //----------------------------------------------Me.------------------------------------------------
                case eOldroydB:

                    switch( m_viscoelastictreatmentType)
                    {
                        case eFullExplicit:

                            m_integrationOps_viscoelasticstress.DefineOdeRhs(&ViscoElasticConstitutive::EvaluateViscoElasticStress_OldroydB_FullExpTreatment_RHS, this);

                            m_integrationOps_viscoelasticstress.DefineImplicitSolve(&ViscoElasticConstitutive::EvaluateViscoElasticStress_OldroydB_FullExpTreatment, this);

                            break;

                        case eImplicitExplicit:

                            m_integrationOps_viscoelasticstress.DefineOdeRhs(&ViscoElasticConstitutive::EvaluateViscoElasticStress_OldroydB_ImpExpTreatment_ExplicitPart, this);

                            m_integrationOps_viscoelasticstress.DefineImplicitSolve(&ViscoElasticConstitutive::EvaluateViscoElasticStress_OldroydB_ImpExpTreatment_Implicitpart, this);

                            break;
                    }

                    break;
            }
        }


    }//End of function void ViscoElasticConstitutive::v_InitObject()

    ViscoElasticConstitutive::~ViscoElasticConstitutive(void)
    {
    }

    void ViscoElasticConstitutive::EvaluateViscoElasticStress_OldroydB_FullExpTreatment(const Array<OneD, const Array<OneD, NekDouble> > &inarray, Array<OneD, Array<OneD, NekDouble> > &outarray, const NekDouble time, const NekDouble aii_Dt)
    {
        //       NekDouble lambda=1.0+aii_Dt*m_ReC/m_We;
        NekDouble lambda=aii_Dt*m_ReC/m_We+1;
        int       physTot = m_fields[0]->GetTotPoints();
        for(int j = 0; j < m_nViscoElasticStressFields; ++j)

        {
            Vmath::Vcopy(physTot,inarray[j],1,outarray[j],1);
            Vmath::Smul(physTot,1/lambda,outarray[j],1,outarray[j],1);
        }
    }

    void ViscoElasticConstitutive::EvaluateViscoElasticStress_OldroydB_FullExpTreatment_RHS(const Array<OneD, const Array<OneD, NekDouble> > &fields, Array<OneD, Array<OneD, NekDouble> > &outarray, const NekDouble time)
    {

        int       i,j,k,nc,coun;
        int       physTot = m_fields[0]->GetTotPoints();

        Array<OneD, NekDouble> wk = Array<OneD, NekDouble>(physTot);
        Array<OneD, NekDouble> tmp = Array<OneD, NekDouble>(physTot);
        Array<TwoD, Array< OneD, NekDouble> > tmp2(m_spacedim,m_spacedim);
        Array<TwoD, Array< OneD, NekDouble> > tmp3(m_spacedim,m_spacedim);
        Array<TwoD, Array< OneD, NekDouble> > tmp4(m_spacedim,m_spacedim);
        Array<OneD, Array< OneD, NekDouble> > F(m_nViscoElasticStressFields);

        for(j=0; j<m_spacedim;j++)
        {
            for(i=0; i<m_spacedim;i++)
            {
                tmp2[i][j]=Array<OneD, NekDouble> (physTot);
                tmp3[i][j]=Array<OneD, NekDouble> (physTot);
                tmp4[i][j]=Array<OneD, NekDouble> (physTot);
            }
        }

        Array<OneD, Array<OneD, NekDouble> >   fields_velocity(m_nConvectiveFields);

        for(i = 0; i < m_nConvectiveFields; ++i)
        {
            fields_velocity[i]  = m_fields[i]->GetPhys();
        }

        // Calculate -V.grad(Tau)
        for(j = 0; j < m_nViscoElasticStressFields; ++j)
        {
            Vmath::Zero(physTot,tmp,1);

            for(i = 0; i < m_nConvectiveFields; ++i)
            {
                m_fields[i]->PhysDeriv(MultiRegions::DirCartesianMap[i],fields[j], wk);
                Vmath::Vvtvp(physTot,wk,1,fields_velocity[i],1,tmp,1,tmp,1);
            }

            Vmath::Smul(physTot,-m_We,tmp,1,outarray[j],1);
        }


        //Calculate grad(V)^T

        for(i = 0; i < m_nConvectiveFields; ++i)
        {
            for(j = 0; j < m_nConvectiveFields; ++j)
            {
                m_fields[i]->PhysDeriv(MultiRegions::DirCartesianMap[i],fields_velocity[j], tmp2[i][j]);
            }
        }

        // put Tau in matrix n by n


        if (m_spacedim == 2){ //2D
            Vmath::Vcopy(physTot,fields[0],1,tmp3[0][0],1);
            Vmath::Vcopy(physTot,fields[1],1,tmp3[0][1],1);
            Vmath::Vcopy(physTot,fields[1],1,tmp3[1][0],1);
            Vmath::Vcopy(physTot,fields[2],1,tmp3[1][1],1);

        }

        else if (m_spacedim == 3){ //3D
            Vmath::Vcopy(physTot,fields[0],1,tmp3[0][0],1);
            Vmath::Vcopy(physTot,fields[1],1,tmp3[0][1],1);
            Vmath::Vcopy(physTot,fields[2],1,tmp3[0][2],1);
            Vmath::Vcopy(physTot,fields[1],1,tmp3[1][0],1);
            Vmath::Vcopy(physTot,fields[3],1,tmp3[1][1],1);
            Vmath::Vcopy(physTot,fields[4],1,tmp3[1][2],1);
            Vmath::Vcopy(physTot,fields[2],1,tmp3[2][0],1);
            Vmath::Vcopy(physTot,fields[4],1,tmp3[2][1],1);
            Vmath::Vcopy(physTot,fields[5],1,tmp3[2][2],1);

        }

        else{
            ASSERTL0(false," not designed for one dimension viscoelastic Stress");
        }


        // Calculate Tau.grad(V)^T

        coun=-1;

        for(i=0; i<m_nConvectiveFields;i++)
        {
            for(j=i; j<m_nConvectiveFields;j++)
            {
                coun=coun+1;

                Vmath::Zero(physTot,wk,1);

                for(k=0; k<m_nConvectiveFields;k++)
                {

                    Vmath::Vvtvp(physTot,tmp3[i][k],1,tmp2[k][j],1,wk,1,wk,1);

                }


                Vmath::Smul(physTot,m_We,wk,1,wk,1);

                Vmath::Vadd(physTot,wk,1,outarray[coun],1,outarray[coun],1);
            }

        }


        //Calculate grad(V)


        for(i = 0; i < m_nConvectiveFields; ++i)

        {

            for(j = 0; j < m_nConvectiveFields; ++j)
            {
                m_fields[i]->PhysDeriv(MultiRegions::DirCartesianMap[j],fields_velocity[i],tmp2[i][j]);

            }
        }



        //Calculate grad(V).Tau


        coun=-1;

        for(i=0; i<m_spacedim;i++)
        {
            for(j=i; j<m_spacedim;j++)
            {
                coun=coun+1;

                Vmath::Zero(physTot,wk,1);

               for(k=0; k<m_spacedim;k++)
                {

                    Vmath::Vvtvp(physTot,tmp2[i][k],1,tmp3[k][j],1,wk,1,wk,1);

                }


                Vmath::Smul(physTot,m_We,wk,1,wk,1);


                Vmath::Vadd(physTot,wk,1,outarray[coun],1,outarray[coun],1);

            }
        }



        //Calculate 2*m_ReM4*D   D=1/2*(grad(V)+grad(V)^T)

        coun=-1;

        for(i = 0; i < m_nConvectiveFields; ++i)

        {
            for(j =i; j < m_nConvectiveFields; ++j)
            {


                coun=coun+1;
                m_fields[i]->PhysDeriv(MultiRegions::DirCartesianMap[j],fields_velocity[i],tmp2[i][j]);
                m_fields[i]->PhysDeriv(MultiRegions::DirCartesianMap[i],fields_velocity[j],tmp3[i][j]);
                Vmath::Vadd(physTot,tmp2[i][j],1,tmp3[i][j],1,wk,1);



                Vmath::Smul(physTot,m_ReM4,wk,1,wk,1);
                Vmath::Vadd(physTot,wk,1,outarray[coun],1,outarray[coun],1);
                Vmath::Smul(physTot,m_ReC/m_We,outarray[coun],1,outarray[coun],1);


            }

        }

    }

    void  ViscoElasticConstitutive::SetUpViscoElasticForcing(Array<OneD, Array<OneD, NekDouble> > &Forcing)
    {


        int i,j;
        int       physTot = m_fields[0]->GetTotPoints();
        int Nx,Ny,n,offset;
        Array<TwoD, Array< OneD, NekDouble> > tmp(m_spacedim,m_spacedim);

        for(j=0; j<m_spacedim;j++)
        {
            for(i=0; i<m_spacedim;i++)
            {
                tmp[i][j]=Array<OneD, NekDouble> (physTot);
            }
        }


        Array<OneD, Array<OneD, NekDouble> >   fields_viscoelasticstress(m_nViscoElasticStressFields);


        for(int n = 0; n < m_nViscoElasticStressFields; ++n)
        {
            fields_viscoelasticstress[n]  = m_fields[n+m_spacedim]->GetPhys();
        }



        if (m_spacedim == 2){ //2D
            Vmath::Vcopy(physTot,fields_viscoelasticstress[0],1,tmp[0][0],1);
            Vmath::Vcopy(physTot,fields_viscoelasticstress[1],1,tmp[0][1],1);
            Vmath::Vcopy(physTot,fields_viscoelasticstress[1],1,tmp[1][0],1);
            Vmath::Vcopy(physTot,fields_viscoelasticstress[2],1,tmp[1][1],1);

        }

        else if (m_spacedim == 3){ //3D


            Vmath::Vcopy(physTot,fields_viscoelasticstress[0],1,tmp[0][0],1);
            Vmath::Vcopy(physTot,fields_viscoelasticstress[1],1,tmp[0][1],1);
            Vmath::Vcopy(physTot,fields_viscoelasticstress[2],1,tmp[0][2],1);
            Vmath::Vcopy(physTot,fields_viscoelasticstress[1],1,tmp[1][0],1);
            Vmath::Vcopy(physTot,fields_viscoelasticstress[3],1,tmp[1][1],1);
            Vmath::Vcopy(physTot,fields_viscoelasticstress[4],1,tmp[1][2],1);
            Vmath::Vcopy(physTot,fields_viscoelasticstress[2],1,tmp[2][0],1);
            Vmath::Vcopy(physTot,fields_viscoelasticstress[4],1,tmp[2][1],1);
            Vmath::Vcopy(physTot,fields_viscoelasticstress[5],1,tmp[2][2],1);

        }

        else{
            ASSERTL0(false," not designed for one dimension viscoelastic Stress");
        }


        for(i = 0; i < m_nConvectiveFields; ++i)
        {

            Vmath::Zero(physTot,Forcing[i],1);
            for(j = 0; j < m_nConvectiveFields; ++j)
            {


                m_fields[i]->PhysDeriv(MultiRegions::DirCartesianMap[j],tmp[i][j], tmp[i][j]);
                Vmath::Vadd(physTot,tmp[i][j],1,Forcing[i],1,Forcing[i],1);

            }


            Vmath::Smul(physTot,m_ReM3,Forcing[i],1,Forcing[i],1);


        }

    }


    void  ViscoElasticConstitutive::SetUpViscoElasticForcing_Weak(Array<OneD, Array<OneD, NekDouble> > &F)
    {


        Array<OneD, Array<OneD, NekDouble> >   fields_viscoelasticstress(m_nViscoElasticStressFields);

        for(int n = 0; n < m_nViscoElasticStressFields; ++n)
        {
            fields_viscoelasticstress[n]  = m_fields[n+m_spacedim]->GetPhys();
        }


        int j;
        int i;

        int  nCoeffs = m_fields[0]->GetNcoeffs();
        int       physTot = m_fields[0]->GetTotPoints();
        Array<OneD, NekDouble> iprod(nCoeffs);

        int ncoeffs = m_fields[0]->GetNcoeffs();

        Array<OneD, Array< OneD, NekDouble> > outarray_weak(m_nConvectiveFields);
        for(j=0;j< m_nConvectiveFields;j++){
            outarray_weak[j]=Array<OneD, NekDouble> (ncoeffs);
        }


        Array<OneD, Array<OneD, NekDouble> >   Forcing(m_nViscoElasticStressFields);
        for(int n = 0; n < m_nConvectiveFields; ++n)
        {
            Forcing[n] = Array<OneD, NekDouble> (ncoeffs);

        }


        Array<TwoD, Array< OneD, NekDouble> > tmp(m_spacedim,m_spacedim);

        for(j=0; j<m_spacedim;j++)
        {
            for(i=0; i<m_spacedim;i++)
            {
                tmp[i][j]=Array<OneD, NekDouble> (physTot);
            }
        }



        if (m_spacedim == 2){ //2D
            Vmath::Vcopy(physTot,fields_viscoelasticstress[0],1,tmp[0][0],1);
            Vmath::Vcopy(physTot,fields_viscoelasticstress[1],1,tmp[0][1],1);
            Vmath::Vcopy(physTot,fields_viscoelasticstress[1],1,tmp[1][0],1);
            Vmath::Vcopy(physTot,fields_viscoelasticstress[2],1,tmp[1][1],1);

        }

        else if (m_spacedim == 3){ //3D

            Vmath::Vcopy(physTot,fields_viscoelasticstress[0],1,tmp[0][0],1);
            Vmath::Vcopy(physTot,fields_viscoelasticstress[1],1,tmp[0][1],1);
            Vmath::Vcopy(physTot,fields_viscoelasticstress[2],1,tmp[0][2],1);
            Vmath::Vcopy(physTot,fields_viscoelasticstress[1],1,tmp[1][0],1);
            Vmath::Vcopy(physTot,fields_viscoelasticstress[3],1,tmp[1][1],1);
            Vmath::Vcopy(physTot,fields_viscoelasticstress[4],1,tmp[1][2],1);
            Vmath::Vcopy(physTot,fields_viscoelasticstress[2],1,tmp[2][0],1);
            Vmath::Vcopy(physTot,fields_viscoelasticstress[4],1,tmp[2][1],1);
            Vmath::Vcopy(physTot,fields_viscoelasticstress[5],1,tmp[2][2],1);

        }

        else{
            ASSERTL0(false," not designed for one dimension viscoelastic Stress");
        }

        for(i = 0; i < m_nConvectiveFields; ++i)
        {

            Vmath::Zero(nCoeffs,F[i],1);
            for(j = 0; j < m_nConvectiveFields; ++j)
            {

                cout<<" test"<< j<<endl;

                m_fields[i]->IProductWRTDerivBase(j, tmp[i][j], iprod);

                cout<<"max F0"<<  Vmath::Vmax(m_fields[0]->GetNcoeffs(),Forcing[0],1)<<endl;
                cout<<"max F1"<<  Vmath::Vmax(m_fields[0]->GetNcoeffs(),Forcing[1],1)<<endl;

                Vmath::Vadd(nCoeffs,iprod,1,F[i], 1, F[i], 1);


            }

            Vmath::Neg( nCoeffs,F[i],1);
            Vmath::Smul(nCoeffs,m_ReM3,F[i],1,F[i],1);

        }
        ViscoElasticConstitutive::AddStressTimesNormalToVelocityNeumannBC(fields_viscoelasticstress,F);

    }

    //========================================== Explicit Implicit Treatment==================================
    void ViscoElasticConstitutive::EvaluateViscoElasticStress_OldroydB_ImpExpTreatment_ExplicitPart(const Array<OneD, const Array<OneD, NekDouble> > &fields, Array<OneD, Array<OneD, NekDouble> > &outarray,  const NekDouble time)
    {

        int i,j;
        int physTot = m_fields[0]->GetTotPoints();
        int ncoeffs = m_fields[0]->GetNcoeffs();
        int Nx,Ny,n,offset;
        Array<OneD, NekDouble> wk = Array<OneD, NekDouble>(physTot);
        Array<OneD, NekDouble> tmp = Array<OneD, NekDouble>(physTot);
        Array<OneD, Array<OneD, NekDouble> >   fields_velocity(m_nConvectiveFields);
        Array<OneD, Array< OneD, NekDouble> > outarray_weak(m_nViscoElasticStressFields);


        for(i = 0; i < m_nConvectiveFields; ++i)
        {
            fields_velocity[i]  = m_fields[i]->GetPhys();
        }

        for(j=0;j<m_nViscoElasticStressFields;j++){
            outarray_weak[j]=Array<OneD, NekDouble> (ncoeffs);
        }


        // Calculate -m_ReC*V.grad(Tau)
        for(j = 0; j < m_nViscoElasticStressFields; ++j)
        {
            Vmath::Zero(physTot,tmp,1);


            for(i = 0; i < m_nConvectiveFields; ++i)
            {

                m_fields[i]->PhysDeriv(MultiRegions::DirCartesianMap[i],fields[j], wk);

                Vmath::Vvtvp(physTot,wk,1,fields_velocity[i],1,tmp,1,tmp,1);

            }


            Vmath::Smul(physTot,-m_ReC,tmp,1,outarray[j],1);

        }


    }// 

    void ViscoElasticConstitutive::EvaluateViscoElasticStress_OldroydB_ImpExpTreatment_Implicitpart(const Array<OneD, const Array<OneD, NekDouble> > &inarray, 
            Array<OneD, Array<OneD, NekDouble> > &outarray, 
            const NekDouble time,
            const NekDouble aii_Dt)
    {


        int i,j;

        int physTot = m_fields[0]->GetTotPoints();
        int ncoeffs = m_fields[0]->GetNcoeffs();
        int Nx,Ny;
        int n,offset;

        Array<TwoD, Array< OneD, NekDouble> > tmp(inarray.num_elements(),inarray.num_elements());
        Array<TwoD, Array< OneD, NekDouble> > gradv(m_spacedim,m_spacedim);
        Array<TwoD, Array< OneD, NekDouble> > G_weak(m_spacedim,m_spacedim);
        Array<OneD, Array< OneD, NekDouble> > rhs(inarray.num_elements());
        Array<OneD, Array< OneD, NekDouble> > rhs_weak(inarray.num_elements());
        Array<OneD, Array< OneD, NekDouble> > outarray_weak(m_nViscoElasticStressFields);



        for(j=0; j<inarray.num_elements();j++)
        {
            rhs[j]=Array<OneD, NekDouble> (physTot);
            rhs_weak[j]=Array<OneD, NekDouble> (ncoeffs);
            for(i=0; i<inarray.num_elements();i++)
            {
                tmp[i][j]=Array<OneD, NekDouble> (physTot);
            }
        }

        for(j=0; j<m_spacedim;j++)
        {
            for(i=0; i<m_spacedim;i++)
            {
                gradv[i][j]=Array<OneD, NekDouble> (physTot);
                G_weak[i][j]=Array<OneD, NekDouble> (ncoeffs);
            }
        }

        for(j=0;j<m_nViscoElasticStressFields;j++){
            outarray_weak[j]=Array<OneD, NekDouble> (ncoeffs);
        }

        Array<OneD, NekDouble> exactsoln(physTot);
        // EvaluateExactSolution(0, exactsoln, m_time);




        Array<OneD, Array<OneD, NekDouble> >   fields_velocity(m_nConvectiveFields);
        for(i = 0; i < m_nConvectiveFields; ++i)
        {
            fields_velocity[i]=Array<OneD, NekDouble> (physTot);
            fields_velocity[i]  = m_fields[i]->GetPhys();
        }


        //Calculate grad(V)


        for(i = 0; i < m_nConvectiveFields; ++i)

        {

            for(j = 0; j < m_nConvectiveFields; ++j)
            {
                m_fields[i]->PhysDeriv(MultiRegions::DirCartesianMap[j],fields_velocity[i],gradv[i][j]);


            }
        }

        //calculate m : the coefficient matrix and m^-1 the inverse matrix

        if (m_spacedim == 2){ //2D
            //m_0,0
            Vmath::Fill(physTot,aii_Dt+m_We/m_ReC,tmp[0][0],1);
            Vmath::Svtvp(physTot,-2*m_We*aii_Dt, gradv[0][0],1,tmp[0][0],1,tmp[0][0],1);



            //m_0,1
            Vmath::Zero(physTot,tmp[0][1],1);
            Vmath::Svtvp(physTot,-2*m_We*aii_Dt, gradv[0][1],1,tmp[0][1],1,tmp[0][1],1);

            //m_0,2
            Vmath::Zero(physTot,tmp[0][2],1);

            //m_1,0
            Vmath::Zero(physTot,tmp[1][0],1);
            Vmath::Svtvp(physTot,-m_We*aii_Dt, gradv[1][0],1,tmp[1][0],1,tmp[1][0],1);

            //m_1,1
            Vmath::Fill(physTot,aii_Dt+m_We/m_ReC,tmp[1][1],1);
            Vmath::Svtvp(physTot,-m_We*aii_Dt, gradv[0][0],1,tmp[1][1],1,tmp[1][1],1);
            Vmath::Svtvp(physTot,-m_We*aii_Dt, gradv[1][1],1,tmp[1][1],1,tmp[1][1],1);


            //m_1,2
            Vmath::Zero(physTot,tmp[1][2],1);
            Vmath::Svtvp(physTot,-m_We*aii_Dt, gradv[0][1],1,tmp[1][2],1,tmp[1][2],1);

            //m_2,0
            Vmath::Zero(physTot,tmp[2][0],1);

            //m_2,1
            Vmath::Zero(physTot,tmp[2][1],1);
            Vmath::Svtvp(physTot,-2*m_We*aii_Dt, gradv[1][0],1,tmp[2][1],1,tmp[2][1],1);

            //m_2,2
            Vmath::Fill(physTot,aii_Dt+m_We/m_ReC,tmp[2][2],1);
            Vmath::Svtvp(physTot,-2.0*m_We*aii_Dt, gradv[1][1],1,tmp[2][2],1,tmp[2][2],1);


        }

        else if (m_spacedim == 3){ //3D


            // m_0,0
            Vmath::Fill(physTot,aii_Dt+m_We/m_ReC,tmp[0][0],1);
            Vmath::Svtvp(physTot,-2*m_We*aii_Dt, gradv[0][0],1,tmp[0][0],1,tmp[0][0],1);

            // m_0,1
            Vmath::Zero(physTot,tmp[0][1],1);
            Vmath::Svtvp(physTot,-2*m_We*aii_Dt, gradv[0][1],1,tmp[0][1],1,tmp[0][1],1);

            // m_0,2
            Vmath::Zero(physTot,tmp[0][2],1);
            Vmath::Svtvp(physTot,-2*m_We*aii_Dt, gradv[0][2],1,tmp[0][2],1,tmp[0][2],1);

            // m_0,3
            Vmath::Zero(physTot,tmp[0][3],1);

            // m_0,4
            Vmath::Zero(physTot,tmp[0][4],1);

            // m_0,5
            Vmath::Zero(physTot,tmp[0][5],1);

            // m_1,0
            Vmath::Zero(physTot,tmp[1][0],1);
            Vmath::Svtvp(physTot,-m_We*aii_Dt, gradv[1][0],1,tmp[1][0],1,tmp[1][0],1);

            //m_1,1
            Vmath::Fill(physTot,aii_Dt+m_We/m_ReC,tmp[1][1],1);
            Vmath::Svtvp(physTot,-m_We*aii_Dt, gradv[0][0],1,tmp[1][1],1,tmp[1][1],1);
            Vmath::Svtvp(physTot,-m_We*aii_Dt, gradv[1][1],1,tmp[1][1],1,tmp[1][1],1);

            // m_1,2
            Vmath::Zero(physTot,tmp[1][2],1);
            Vmath::Svtvp(physTot,-m_We*aii_Dt, gradv[1][2],1,tmp[1][2],1,tmp[1][2],1);

            // m_1,3
            Vmath::Zero(physTot,tmp[1][3],1);
            Vmath::Svtvp(physTot,-m_We*aii_Dt, gradv[0][1],1,tmp[1][3],1,tmp[1][3],1);

            // m_1,4
            Vmath::Zero(physTot,tmp[1][4],1);
            Vmath::Svtvp(physTot,-m_We*aii_Dt, gradv[0][2],1,tmp[1][4],1,tmp[1][4],1);

            // m_1,5
            Vmath::Zero(physTot,tmp[1][5],1);

            // m_2,0
            Vmath::Zero(physTot,tmp[2][0],1);
            Vmath::Svtvp(physTot,-m_We*aii_Dt, gradv[2][0],1,tmp[2][0],1,tmp[2][0],1);


            // m_2,1
            Vmath::Zero(physTot,tmp[2][1],1);
            Vmath::Svtvp(physTot,-m_We*aii_Dt, gradv[2][1],1,tmp[2][1],1,tmp[2][1],1);

            //m_2,2
            Vmath::Fill(physTot,aii_Dt+m_We/m_ReC,tmp[2][2],1);
            Vmath::Svtvp(physTot,-m_We*aii_Dt, gradv[0][0],1,tmp[2][2],1,tmp[2][2],1);
            Vmath::Svtvp(physTot,-m_We*aii_Dt, gradv[2][2],1,tmp[2][2],1,tmp[2][2],1);

            // m_2,3
            Vmath::Zero(physTot,tmp[2][3],1);

            // m_2,4
            Vmath::Zero(physTot,tmp[2][4],1);
            Vmath::Svtvp(physTot,-m_We*aii_Dt, gradv[0][1],1,tmp[2][4],1,tmp[2][4],1);

            // m_2,5
            Vmath::Zero(physTot,tmp[2][5],1);
            Vmath::Svtvp(physTot,-m_We*aii_Dt, gradv[0][2],1,tmp[2][5],1,tmp[2][5],1);

            // m_3,0
            Vmath::Zero(physTot,tmp[3][0],1);

            //m_3,1
            Vmath::Zero(physTot,tmp[3][1],1);
            Vmath::Svtvp(physTot,-2*m_We*aii_Dt, gradv[1][0],1,tmp[3][1],1,tmp[3][1],1);

            // m_3,2
            Vmath::Zero(physTot,tmp[3][2],1);


            //m_3,3
            Vmath::Fill(physTot,aii_Dt+m_We/m_ReC,tmp[3][3],1);
            Vmath::Svtvp(physTot,-2*m_We*aii_Dt, gradv[1][1],1,tmp[3][3],1,tmp[3][3],1);


            //m_3,4
            Vmath::Zero(physTot,tmp[3][4],1);
            Vmath::Svtvp(physTot,-2*m_We*aii_Dt, gradv[1][2],1,tmp[3][4],1,tmp[3][4],1);

            // m_3,5
            Vmath::Zero(physTot,tmp[3][5],1);

            // m_4,0
            Vmath::Zero(physTot,tmp[4][0],1);

            //m_4,1
            Vmath::Zero(physTot,tmp[4][1],1);
            Vmath::Svtvp(physTot,-m_We*aii_Dt, gradv[2][0],1,tmp[4][1],1,tmp[4][1],1);

            //m_4,2
            Vmath::Zero(physTot,tmp[4][2],1);
            Vmath::Svtvp(physTot,-m_We*aii_Dt, gradv[1][0],1,tmp[4][2],1,tmp[4][2],1);

            //m_4,3
            Vmath::Zero(physTot,tmp[4][3],1);
            Vmath::Svtvp(physTot,-m_We*aii_Dt, gradv[2][1],1,tmp[4][3],1,tmp[4][3],1);

            //m_4,4
            Vmath::Fill(physTot,aii_Dt+m_We/m_ReC,tmp[4][4],1);
            Vmath::Svtvp(physTot,-m_We*aii_Dt, gradv[1][1],1,tmp[4][4],1,tmp[4][4],1);
            Vmath::Svtvp(physTot,-m_We*aii_Dt, gradv[2][2],1,tmp[4][4],1,tmp[4][4],1);

            //m_4,5
            Vmath::Zero(physTot,tmp[4][5],1);
            Vmath::Svtvp(physTot,-m_We*aii_Dt, gradv[1][2],1,tmp[4][5],1,tmp[4][5],1);

            // m_5,0
            Vmath::Zero(physTot,tmp[5][0],1);

            // m_5,1
            Vmath::Zero(physTot,tmp[5][1],1);

            //m_5,2
            Vmath::Zero(physTot,tmp[5][2],1);
            Vmath::Svtvp(physTot,-2*m_We*aii_Dt, gradv[2][0],1,tmp[5][2],1,tmp[5][2],1);

            // m_5,3
            Vmath::Zero(physTot,tmp[5][3],1);

            //m_5,4
            Vmath::Zero(physTot,tmp[5][4],1);
            Vmath::Svtvp(physTot,-2*m_We*aii_Dt, gradv[2][1],1,tmp[5][4],1,tmp[5][4],1);

            //m_5,5
            Vmath::Fill(physTot,aii_Dt+m_We/m_ReC,tmp[5][5],1);
            Vmath::Svtvp(physTot,-2*m_We*aii_Dt, gradv[2][2],1,tmp[5][5],1,tmp[5][5],1);



        }

        else{

            ASSERTL0(false," not designed for one dimension viscoelastic Stress");
        }


        //callculate rhs: (1-Rmu)(gardv+gradv^T)+m_We/m_ReC*tau_intermediate


        int      coun=-1;

        for(i = 0; i < m_nConvectiveFields; ++i)

        {
            for(j =i; j < m_nConvectiveFields; ++j)
            {
                coun=coun+1;
                Vmath::Zero(physTot,rhs[coun],1);
                Vmath::Zero(ncoeffs,rhs_weak[coun],1);

                m_fields[i]->PhysDeriv(MultiRegions::DirCartesianMap[j],fields_velocity[i],gradv[i][j]);


                m_fields[i]->PhysDeriv(MultiRegions::DirCartesianMap[i],fields_velocity[j],rhs[coun]);


                Vmath::Vadd(physTot,rhs[coun],1,gradv[i][j],1,rhs[coun],1);
                Vmath::Smul(physTot,aii_Dt*m_ReM4,rhs[coun],1,rhs[coun],1);

            }

        }


        coun=-1;

        for(i = 0; i < m_nConvectiveFields; ++i)

        {
            for(j =i; j < m_nConvectiveFields; ++j)
            {
                coun=coun+1;
                Vmath::Svtvp(physTot,m_We/m_ReC, inarray[coun],1,rhs[coun],1,rhs[coun],1);
            }
        }


        // tau=(m^-1).rhs


        if (m_spacedim == 2){ //2D


            for(i=0; i<physTot;i++)
            {

                double buf1[]= {tmp[0][0][i],tmp[1][0][i],tmp[2][0][i],
                    tmp[0][1][i],tmp[1][1][i],tmp[2][1][i],
                    tmp[0][2][i],tmp[1][2][i],tmp[2][2][i]};


                double buf2[]= {rhs[0][i],rhs[1][i],rhs[2][i]};

                NekMatrix<double> m(3,3,buf1);

                m.Invert();

                NekVector<double> RHS(3, buf2);
                NekVector<double> result = m*RHS;

                for(j=0; j<3;j++)
                    outarray[j][i]=result(j);



            }
        }

        else if (m_spacedim == 3){ //3D

            for(i=0; i<physTot;i++)
            {

                double buf1[]= {tmp[0][0][i],tmp[1][0][i],tmp[2][0][i],tmp[3][0][i],tmp[4][0][i],tmp[5][0][i],
                    tmp[0][1][i],tmp[1][1][i],tmp[2][1][i],tmp[3][1][i],tmp[4][1][i],tmp[5][1][i],
                    tmp[0][2][i],tmp[1][2][i],tmp[2][2][i],tmp[3][2][i],tmp[4][2][i],tmp[5][2][i],
                    tmp[0][3][i],tmp[1][3][i],tmp[2][3][i],tmp[3][3][i],tmp[4][3][i],tmp[5][3][i],
                    tmp[0][4][i],tmp[1][4][i],tmp[2][4][i],tmp[3][4][i],tmp[4][4][i],tmp[5][4][i],
                    tmp[0][5][i],tmp[1][5][i],tmp[2][5][i],tmp[3][5][i],tmp[4][5][i],tmp[5][5][i]};

                double buf2[]= {rhs[0][i],rhs[1][i],rhs[2][i],rhs[3][i],rhs[4][i],rhs[5][i]};

                NekMatrix<double> m(6,6,buf1);
                m.Invert();

                NekVector<double> RHS(6, buf2);
                NekVector<double> result = m*RHS;

                for(j=0; j<6;j++)
                    outarray[j][i]=result(j);

            }

        }

        else{
            ASSERTL0(false," not designed for one dimension viscoelastic Stress");
        }

        for(j = 0; j < m_nViscoElasticStressFields; ++j)
        {

            m_fields[j+m_spacedim]->FwdTrans(outarray[j],outarray_weak[j]);
            m_fields[j+m_spacedim]->BwdTrans(outarray_weak[j],outarray[j]);
        }



    }//
    //-------------------------------Me.------------------------------

    void ViscoElasticConstitutive::EvaluateAggregateSize_FullExpTreatment_RHS(const Array<OneD, const Array<OneD, NekDouble> > &fields, Array<OneD, Array<OneD, NekDouble> > &outarray, const NekDouble time)
    {
        // ........variables........ 
        int i, j;
        int       physTot = m_fields[0]->GetTotPoints();
        Array<TwoD, Array< OneD, NekDouble> > gradv(m_spacedim,m_spacedim);
        Array<TwoD, Array< OneD, NekDouble> > gradvT(m_spacedim,m_spacedim);
        Array<TwoD, Array< OneD, NekDouble> > Gammadot(m_spacedim,m_spacedim);
        Array<OneD, NekDouble> gammadot = Array<OneD, NekDouble> (physTot);
        Array<OneD, Array<OneD, NekDouble> >   fields_velocity(m_nConvectiveFields);

        for(i = 0; i < m_nConvectiveFields; ++i)
        {
            fields_velocity[i]  = m_fields[i]->GetPhys();
        }
        //*/ is these lines of code necessary? (ask about)
        for(j=0; j<m_spacedim;j++)
        {
            for(i=0; i<m_spacedim;i++)
            {
                gradv[i][j] = Array<OneD, NekDouble> (physTot);
                gradvT[j][i] = Array<OneD, NekDouble> (physTot);    // gradient of v transpose
                Gammadot[i][j] = Array<OneD, NekDouble> (physTot);
            }
        }
        //*/        
        double gammadotCr, gammadotMax;
        Array<OneD, NekDouble> a = Array<OneD, NekDouble> (physTot);
        Array<OneD, NekDouble> b = Array<OneD, NekDouble> (physTot);
        double a13, a12, a10, a23, a22, a21, a20;
        double teta, beta, etaInf, etaZero;
        Array<OneD, NekDouble> nst = Array<OneD, NekDouble> (physTot);
        Array<OneD, NekDouble> AggRhs = Array<OneD, NekDouble> (physTot);
        int m; // maybe it's not supposed to be an integer (ask about)
        double N0, landaH;
        //        double gn, muNew;
       
        //calculate Gammadot type: 2D matrix of NekDouble-----------------------------------------------------
        for(i = 0; i < m_nConvectiveFields; ++i)
        {
            for(j = 0; j < m_nConvectiveFields; ++j)
            {
                m_fields[i]->PhysDeriv(MultiRegions::DirCartesianMap[j],fields_velocity[i],gradv[i][j]);
                m_fields[i]->PhysDeriv(MultiRegions::DirCartesianMap[i],fields_velocity[j],gradvT[i][j]);
                Vmath::Vadd(physTot,gradv[i][j],1,gradvT[i][j],1,Gammadot[i][j],1);
                Vmath::Smul(physTot,0.5,Gammadot[i][j],1,Gammadot[i][j],1); // is this line necessary?

            }
        }
       //calculate gammadot type:1D vector of Nekdouble-------------------------------------------------------
        Vmath::Zero(physTot, gammadot,1);
        for(i = 0; i < m_nConvectiveFields; ++i)
        {
            for(j =0; j < m_nConvectiveFields; ++j)
            {
                Vmath::Vvtvp(physTot,Gammadot[i][j],1,Gammadot[i][j],1,gammadot,1,gammadot,1);
            }
        }

        // Calculate shear rate: \dot{\gamma}= \sqrt[+]{\frac{1}{2}\dot\gamma:\dot\gamma}

        Vmath::Smul(physTot,1.0/2.0,gammadot,1,gammadot,1);

        for(i=0; i<physTot;i++)
        {
            gammadot[i]=sqrt(gammadot[i]);
        }

        //
        //calculate a
        Vmath::Zero(physTot, a,1);

        for(i=0; i<physTot;i++)
        {
            if(gammadot[i] <= gammadotCr )
            {
                a[i] = a13*pow(gammadot[i],3)+a12*pow(gammadot[i],2)+a10;
            }
            else if(gammadot[i] <= gammadotMax)
            {
                a[i] = a23*pow(gammadot[i],3)+a22*pow(gammadot[i],2)+a21*gammadot[i]+a20;
            }
        }
        //calculate teta
        teta = beta*etaInf/etaZero;
        //calculate nst
        //edit by pow in Vmath
        for(i=0; i<physTot;i++)
        {
            nst[i] = etaZero/etaInf*(1+teta*pow(gammadot[i],m))/(1+beta*pow(gammadot[i],m))*(1+3/2*a[i]*N0*landaH); 
        }
        //calculate b
        for(i=0; i<physTot;i++)
        {
            b[i] = a[i]*N0/nst[i]/(nst[i]-1);
        }
        //evaluate n
        for(i=0; i<physTot;i++)
        {
            // combine these two lines
            AggRhs[i] = 1.0/2.0*b[i]*(nst[i]-fields[0][i])*(fields[0][i]+nst[i]-1);
        }

        Vmath::Vcopy(physTot, AggRhs, 1, outarray[0], 1);

        //calculate gn
        //        gn = 0.5*b*(n-1)+a*N0/n;
        //calcullate mu and do something for the mu update
        //        muNew = landaH*n/(1+gn*n*landaH);
        printf("function1\n");
        //        exit(-1);
    }

    void ViscoElasticConstitutive::EvaluateAggregateSize_FullExpTreatment(const Array<OneD, const Array<OneD, NekDouble> > &inarray, Array<OneD, Array<OneD, NekDouble> > &outarray, const NekDouble time, const NekDouble aii_Dt)
    {
        int       physTot = m_fields[0]->GetTotPoints();
        //*/
         
            Vmath::Vcopy(physTot,inarray[0],1,outarray[0],1);

        //*/
        printf("function2\n");
        //        exit(-1);
    }
    //-----------------------------------------------------Me.--------------------------------------------------

    // =====================================Weak formulation==================================================
    // deleted from here


    void  ViscoElasticConstitutive::AddStressTimesNormalToVelocityNeumannBC(const Array<OneD, const Array<OneD, NekDouble> > &stressin, Array<OneD, Array<OneD, NekDouble> > &Weakoutarray)
    {

        int i,el,n,cnt, offset, phys_offset,nq;
        int edge, elmtid;
        Array<OneD, NekDouble> e_outarray;
        Array<OneD, const NekDouble> txx,txy,tyy,P;
        // Array<OneD, Array<OneD, StdRegions::StdExpansion1DSharedPtr> >
        //    elmtToTrace = m_traceMap->GetElmtToTrace();
        StdRegions::StdExpansion1DSharedPtr EdgeExp;

        Array<OneD, int> ElmtID,EdgeID;

        //Fills ElmtID and EdgeID with global ids of elements id number and edges id number of the boundary expansion
        m_fields[0]->GetBoundaryToElmtMap(ElmtID,EdgeID);

        int maxpts = 0;
        int ncoeffs = m_fields[0]->GetNcoeffs();
        int phystot = m_fields[0]->GetTotPoints();
        Array<OneD, NekDouble> pressurecoeff(ncoeffs,0.0);
        Array<OneD, NekDouble> pressurephys(phystot,0.0);

        // Project pressure onto velocity space
        m_pressure->BwdTrans(m_pressure->GetCoeffs(), m_pressure->UpdatePhys());
        // ask about the line below ( without enforcing BCs)
        // project onto velocity space without enforcing BCs
        m_fields[0]->FwdTrans_IterPerExp(m_pressure->GetPhys(),pressurecoeff);
        m_fields[0]->BwdTrans(pressurecoeff,pressurephys);


        // find the maximum values of points
        for(cnt = n = 0; n < m_fields[0]->GetBndConditions().num_elements(); ++n)
        {
            //Loop over all elements in boundary region
            //cnt starts with element 0 then going through elements
            for(i = 0; i < m_fields[0]->GetBndCondExpansions()[n]->GetExpSize(); ++i)
            {
                maxpts = max(maxpts, m_fields[0]->GetExp(ElmtID[cnt++])->GetTotPoints());
            }
        }

        Array<OneD, NekDouble> txxedge(3*maxpts);
        Array<OneD, NekDouble> txyedge = txxedge + maxpts;
        Array<OneD, NekDouble> tyyedge = txyedge + maxpts;

        // Neumann BC for velocity component u
        for(cnt = n = 0; n < m_fields[0]->GetBndConditions().num_elements(); ++n)
        {

            // Waters solution
            if((m_fields[0]->GetBndConditions()[n])->GetBoundaryConditionType() == SpatialDomains::eNeumann)
            {

                // loop over elements along boundary
                for(el = 0; el < m_fields[0]->GetBndCondExpansions()[n]->GetExpSize(); ++el,cnt++)
                {



                    EdgeExp =  boost::dynamic_pointer_cast<StdRegions::StdExpansion1D> (m_fields[0]->GetBndCondExpansions()[n]->GetExp(el));

                    nq   = m_fields[0]->GetExp(elmtid)->GetTotPoints();

                    // grab edge on the boundary
                    edge = EdgeID[cnt];
                    elmtid = ElmtID[cnt];

                    offset = m_fields[0]->GetPhys_Offset(elmtid);

                    txx = stressin[0] + offset;
                    txy = stressin[1] + offset;
                    //P = (m_pressure->GetPhys())+ m_pressure->GetPhys_Offset(elmtid);
                    P = pressurephys + offset;

                    Array<OneD, NekDouble> cauchyxx(nq,0.0);

                    Vmath::Vsub(nq,txx,1,P,1,cauchyxx,1);
                    m_fields[0]->GetExp(elmtid)->GetEdgePhysVals(edge,EdgeExp,cauchyxx,txxedge);
                    m_fields[0]->GetExp(elmtid)->GetEdgePhysVals(edge,EdgeExp,txy,txyedge);

                    int nquad_e = EdgeExp->GetNumPoints(0);

                    e_outarray = Weakoutarray[0] + m_fields[0]->GetCoeff_Offset(elmtid);

                }
            }


            cnt +=m_fields[0]->GetBndCondExpansions()[n]->GetExpSize();

        }



        // Neumann BC for velocity component u
        for(cnt = n = 0; n < m_fields[0]->GetBndConditions().num_elements(); ++n)
        {
            // Waters solution
            if((m_fields[0]->GetBndConditions()[n])->GetBoundaryConditionType() == SpatialDomains::eNeumann)
            {


                // loop over elements along boundary
                for(el = 0; el < m_fields[0]->GetBndCondExpansions()[n]->GetExpSize(); ++el,cnt++)
                {
                    EdgeExp =  boost::dynamic_pointer_cast<StdRegions::StdExpansion1D> (m_fields[0]->GetBndCondExpansions()[n]->GetExp(el));

                    nq   = m_fields[0]->GetExp(elmtid)->GetTotPoints();
                    // grab edge on the boundary
                    edge = EdgeID[cnt];
                    elmtid = ElmtID[cnt];

                    offset = m_fields[0]->GetPhys_Offset(elmtid);

                    txy = stressin[1] + offset;
                    tyy = stressin[2] + offset;
                    //P = (m_pressure->GetPhys())+ m_pressure->GetPhys_Offset(elmtid);
                    P = pressurephys + offset;

                    Array<OneD, NekDouble> cauchyyy(nq,0.0);
                    Vmath::Vsub(nq,txx,1,P,1,cauchyyy,1);

                    m_fields[0]->GetExp(elmtid)->GetEdgePhysVals(edge,EdgeExp,txy,txyedge);
                    m_fields[0]->GetExp(elmtid)->GetEdgePhysVals(edge,EdgeExp,cauchyyy,tyyedge);

                    int nquad_e = EdgeExp->GetNumPoints(0);

                    e_outarray = Weakoutarray[1] + m_fields[0]->GetCoeff_Offset(elmtid);


                }
            }

            cnt +=m_fields[0]->GetBndCondExpansions()[n]->GetExpSize();
        }

    }

}// end of name space

//Revision 1.0 2015/06/23  Oldroyd-B model Jafari
