using System;

namespace HaploGrepSharp
{
    public enum Mutations {
        A, G, C, T, INS, DEL
    }
    
    public static class MutationAssigner {        
        public static Mutations getBase(string mutation)
        {

            if ((mutation.Equals("A")) || (mutation.Equals("a")))

                return Mutations.A;

            if ((mutation.Equals("C")) || (mutation.Equals("c")))

                return Mutations.C;

            if ((mutation.Equals("G")) || (mutation.Equals("g")))

                return Mutations.G;

            if ((mutation.Equals("T")) || (mutation.Equals("t")))
            {
                return Mutations.T;
            }
          
            throw new HaploGrepException("This base does not match a valid base: " + mutation);
        }
        public static string getBase(Mutations mutant)
        {
            switch (mutant)
            {
                case Mutations.A:
                    return "A";
                case Mutations.C:
                    return "C";
                case Mutations.G:
                    return "G";
                case Mutations.T:
                    return "T";
                default:
                    throw new HaploGrepException("Tried to get base for non A,C,G,T mutation");
            }
        }
        public static bool MutationIsBasePair(Mutations mutation)
        {
            return mutation == Mutations.A || mutation == Mutations.C || mutation == Mutations.G || mutation == Mutations.T;
        }
        public static bool MutationIsComplex(Mutations mut)
        {
            return mut == Mutations.INS || mut == Mutations.DEL; 
        }
    }
}

