class Solution {
    public String solution(String a, String b) {
        int carry = 0;
        String answer = "";
        int i = 0;
        while(i < Math.max(a.length(), b.length())){
            int num1 = 0 <= a.length() - (i + 1) ? a.charAt(a.length() - (i + 1)) - '0' : 0;
            int num2 = 0 <= b.length() - (i + 1) ? b.charAt(b.length() - (i + 1)) - '0' : 0;
            int sum = num1 + num2 + carry;
            carry = sum / 10;
            answer = sum % 10 + answer;
            i++;
        }
        if(carry > 0) answer = carry + answer;
        return answer;
    }
}