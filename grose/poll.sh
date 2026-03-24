#!/usr/bin/env bash

URL="https://api.binance.com/api/v3/ticker/price?symbol=BTCUSDT"
while true; do
    curl -s "$URL" | jq -r '.price'
    sleep 0.5
done
