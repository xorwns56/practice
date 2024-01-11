import java.util.*;
class Solution {
    public int[] solution(int[] arr) {
        int new_length = 1;
        while(new_length < arr.length) new_length <<= 1;
        return Arrays.copyOfRange(arr, 0, new_length);
    }
}