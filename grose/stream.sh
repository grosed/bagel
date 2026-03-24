#!/usr/bin/env bash

websocat wss://stream.binance.com:9443/ws/btcusdt@trade | jq -r '.p'
