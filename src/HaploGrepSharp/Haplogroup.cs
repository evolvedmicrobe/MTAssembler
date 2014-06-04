using System;
using System.Collections;
using System.Collections.Generic;
using System.Text.RegularExpressions;

namespace HaploGrepSharp
{
    public class Haplogroup
    {
        //TODO: Rename later
        internal List<object> subParts = new List<object>();
        public string id {get;internal set;}
        public static Regex p = new Regex("\\d+|'", RegexOptions.Compiled);

        public Haplogroup(string _haploGroup)
        {
            this.id = _haploGroup;
            this.subParts.Add(_haploGroup);
        }
        public void changeHaplogroupFormat(string haploGroup)
        {
            int position = 0;            
            bool toggle = true;
            int pos = haploGroup.IndexOf("+");
            string specialEnd = "";
            if (pos > -1)
            {
                specialEnd = haploGroup.Substring(pos);
                haploGroup = haploGroup.Substring(0, pos);
            }
            while (position < haploGroup.Length)
            {
                if (toggle)
                {
                    char c = haploGroup[position];
                    if ((haploGroup.Length > position + 1) && (c == 'H') && (haploGroup[position + 1] == 'V'))
                    {
                        this.subParts.Add("HV");
                        position += 2;
                        toggle = false;

                    }
                    else if ((haploGroup.Length > position + 1) && (c == 'C') && (haploGroup[position + 1] == 'Z'))
                    {
                        this.subParts.Add("CZ");
                        position += 2;
                        toggle = false;

                    }
                    else
                    {
                        this.subParts.Add(Convert.ToChar(haploGroup[position]));
                        position++;
                        toggle = false;

                    }

                }
                else
                {
                    throw new Exception("You need to translate this");
                    //p.Match((haploGroup.Substring(position, haploGroup.Length - position)
                    //Matcher m = p.matcher(haploGroup.Substring(position, haploGroup.Length - position));

                    //if (m.find())
                    //{
                    //    char? c = Convert.ToChar(haploGroup[position + m.start()]);
                    //    if (c.ToString().Equals("'"))
                    //    {
                    //        this.subParts.Add("'");

                    //    }
                    //    else
                    //    {
                    //        this.subParts.Add(Convert.ToInt32(Convert.ToInt32(haploGroup.Substring(position + m.start(), position + m.end() - (position + m.start())))));

                    //    }
                    //    position += m.end();
                    //    toggle = true;

                    //}
                    //else
                    //{
                    //    this.subParts.Add("'");
                    //    toggle = true;

                    //}

                }

            }
            if (pos > -1)
            {
                this.subParts.Add(specialEnd);
            }

        }
        public override bool Equals(object haplogroup)
        {
            if (!(haplogroup is Haplogroup))
            {
                return false;
            } 
            Haplogroup c = (Haplogroup)haplogroup;
            if (!this.id.Equals(c.id))
            {
                return false;

            }
            return true;

        }
        public bool isSuperHaplogroup(Haplogroup hgToCheck)
        {
            if (!(hgToCheck is Haplogroup))
            {
                return false;

            } Haplogroup c = hgToCheck;
            if (!c.id.Contains(this.id))
            {
                return false;

            }
            return true;

        }
        public override string ToString()
        {
            return String.Join("", subParts);
        }
        public override int GetHashCode()
        {
            throw new NotImplementedException();
            return base.GetHashCode();
        }

        public string AccessionId;
        public string Reference;
    }

}