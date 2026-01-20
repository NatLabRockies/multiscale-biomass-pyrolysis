ps.updateHTC();


if(ps.dict().lookupOrDefault<Switch>("useOperatorSplittingInEEqn",true))
{
    #include "operatorSplittingEEqn.H"
}
else
{
    #include "implicitCouplingEEqn.H"
}
