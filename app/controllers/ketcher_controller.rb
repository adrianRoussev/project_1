class KetcherController < ApplicationController
    def env_url
      render file: Rails.root.join('public', 'ketcher', 'standalone', 'index.html'), layout: false
    end
end
