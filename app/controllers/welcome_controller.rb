class WelcomeController < ApplicationController

    def index
      @proteins= Protein.all
      @inhibitors=Inhibitor.all
    end

  end
  