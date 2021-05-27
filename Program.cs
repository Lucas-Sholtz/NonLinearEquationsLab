using System;

namespace NonLinearEquationsLab
{
    class Program
    {
        static void Main(string[] args)
        {
            //NLE.RunNewtonMethod(5);
            //NLE.RunNewtonMethod(NLE.GetTwoFunctions());
            //NLE.RunModifiedNewtonMethod(5);
            //NLE.RunModifiedNewtonMethod(NLE.GetTwoFunctions());
            NLE.RunRelaxationMethod(10);
            //NLE.RunRelaxationMethod(NLE.GetTwoFunctions());
        }
    }
}
