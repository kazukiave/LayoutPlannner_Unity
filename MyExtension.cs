using System.Collections;
using System.Collections.Generic;
using UnityEngine;

namespace MyExtension
{
    public  enum Child
    {
        A = 0,
        B = 1,
    }

    public enum Axis
    {
        Up = 0,
        Right,
        Down,
        Left,
        None,
    }

    public enum eStatus
    {
        Alive,
        Dead,
    }

    public static class MyExtension
    {
     

        public static void Jitter<T>(this IList<T> list)
        {
            for (int i = list.Count - 1; i > 0; i--)
            {
                int j = Random.Range(0, i + 1);
                var tmp = list[i];
                list[i] = list[j];
                list[j] = tmp;
            }
        }

        public static Axis ToAxis<T>(this int axisInt)
        {
            if (axisInt == (int)Axis.Right)
            {
                return Axis.Right;
            }

            if (axisInt == (int)Axis.Left)
            {
                return Axis.Left;
            }

            if (axisInt == (int)Axis.Up)
            {
                return Axis.Up;
            }

            if (axisInt == (int)Axis.Down)
            {
                return Axis.Down;
            }

            return Axis.None;

        }
    }
}
