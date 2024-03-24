class Solution {
    public String solution(int n) {
        int[] arr = new int[] { 1, 2, 4 };
        int digit = 1;
        while(n - Math.pow(arr.length, digit) > 0) n -= Math.pow(arr.length, digit++);
        n--;
        StringBuilder stringBuilder = new StringBuilder();
        while(digit > 0){
            stringBuilder.insert(0, arr[n % arr.length]);
            n /= arr.length;
            digit--;
        }
        return stringBuilder.toString();
    }
    
}