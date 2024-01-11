class Solution {
    public int solution(String binomial) {
        String[] sp = binomial.split("\\s");
        int a = Integer.parseInt(sp[0]);
        int b = Integer.parseInt(sp[2]);
        char op = sp[1].charAt(0);
        return op == '+' ? a + b : (op == '-' ? a - b : (a * b));
    }
}