using System;

namespace NonLinearEquationsLab
{
    class Program
    {
        static void Main(string[] args)
        {
            NLE.RunNewtonMethod(NLE.GetNFunctions(10));
            //NLE.RunModifiedNewtonMethod(NLE.GetTwoFunctions());
        }
    }
}
