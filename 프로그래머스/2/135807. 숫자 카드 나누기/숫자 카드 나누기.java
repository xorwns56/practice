class Solution {
    public int solution(int[] arrayA, int[] arrayB) {
        int gcdA = arrayA[0];
        int gcdB = arrayB[0];
        for(int i = 1; i < arrayA.length; i++) gcdA = gcd(gcdA, arrayA[i]);
        for(int i = 1; i < arrayB.length; i++) gcdB = gcd(gcdB, arrayB[i]);
        int answerA = gcdA;
        for(int i = 0; i < arrayB.length; i++){
            if(arrayB[i] % gcdA == 0){
                answerA = 0;
                break;
            }
        }
        int answerB = gcdB;
        for(int i = 0; i < arrayA.length; i++){
            if(arrayA[i] % gcdB == 0){
                answerB = 0;
                break;
            }
        }
        return Math.max(answerA, answerB);
    }
    public int gcd(int a, int b){
        if(a == 0) return b;
        return gcd(b % a, a);
    }
}