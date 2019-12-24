/*
Say you have an array for which the ith element is the price of a given stock on day i.
If you were only permitted to complete at most one transaction (i.e., buy one and sell one share of the stock), design an algorithm to find the maximum profit.
Note that you cannot sell a stock before you buy one.

Example 1:
Input: [7,1,5,3,6,4]
Output: 5
Explanation: Buy on day 2 (price = 1) and sell on day 5 (price = 6), profit = 6-1 = 5.
             Not 7-1 = 6, as selling price needs to be larger than buying price.

Example 2:
Input: [7,6,4,3,1]
Output: 0
Explanation: In this case, no transaction is done, i.e. max profit = 0.
*/

class Solution {
    public int maxProfit(int[] prices) {
        int maxProfit = 0; //최대 이익을 0으로 초기화
        for(int i=0;i<prices.length-1;i++){ //구매한 날짜의 인덱스 i, 마지막 날의 전날까지만 구매하도록 함
            for(int j=i+1;j<prices.length;j++){ //판매하는 날짜의 인덱스 j, 최소한 구매한 날짜의 다음 날에 판매해야 함
                int profit = prices[j]-prices[i];
                if(maxProfit < profit) maxProfit = profit;
            }
        }
        return maxProfit;
    }
}
