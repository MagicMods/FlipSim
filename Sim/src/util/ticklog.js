export class TickLog {
    constructor(timer) {
        this.timer = timer;
        this.tick = false;
        this.lastTick = 0;
    }

    update() {
        const now = Date.now();
        if (now - this.lastTick >= this.timer) {
            this.tick = true;
            this.lastTick = now;
        }

    }

    GetTick() {
        return this.tick;
    }

    ResetTick() {
        this.tick = false;
    }
}

